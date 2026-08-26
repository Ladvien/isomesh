# The terrain array — where these two PNGs came from

`game_dig` textures its terrain with a triplanar blend of the two files in this
directory. They are not the artists' files: they are **four material sets**,
each of them four maps, packed into two RGBA **texture arrays** — because WebGL2
counts bound textures and `StandardMaterial` already spends six of them, and
because sampling four sets from eight separate images would spend eight more.
Two array images cost two bindings no matter how many layers they hold. This
file is the record of what was packed and how, so the pack can be redone or a
different set swapped in without guessing.

The bytes are committed rather than fetched at build time. The Open SciVis
volumes `.gitignore` excludes (M-006) have a public URL to fetch from; these
sets do not, and the repository already tracks comparable binaries — the largest
tracked file is `docs/gifs/flying-through-the-rock.gif` at 7.4 MB.

## The four layers, and why in this order

The layer index **is** the array slice index, and `triplanar.wgsl`'s
`LAYER_GRASS`/`LAYER_DIRT_SURFACE`/`LAYER_DIRT_DEEP` constants and
`game_dig.rs`'s `LAYER_CONCRETE` are that index written down a second time. The
order below is therefore load-bearing: `Image::reinterpret_stacked_2d_as_array`
slices a stacked image **top-down**, so an inverted stack compiles, renders, and
paints the boundary walls with grass.

| layer | what | source archive | role in the shader |
|---|---|---|---|
| 0 | grass | `Grass002_1K-JPG.zip` | up-facing terrain, `n.y` above ~0.82 |
| 1 | surface dirt with leaves | `Ground023_1K-JPG.zip` | slopes and pit lips |
| 2 | deep dirt | `Ground051_2K-JPG.zip` | anything below `y ≈ -1.6`, i.e. inside a tunnel |
| 3 | damaged concrete | `ground_0046_1k_38H6pY.zip` | the five sandbox walls, forced with `settings.z = 3` |

Layers 0–2 are blended per fragment by slope and world height; layer 3 is never
blended into them. It is selected outright by a second material instance, which
is what makes the sandbox boundary read as built rather than grown.

## Sources

All under `/mnt/codex_fs/game_assets/textures/pbr/`.

| archive | bytes | origin |
|---|---|---|
| `Grass002_1K-JPG.zip` | 10,641,927 | ambientCG |
| `Ground023_1K-JPG.zip` | 10,511,615 | ambientCG |
| `Ground051_2K-JPG.zip` | 28,035,367 | ambientCG |
| `ground_0046_1k_38H6pY.zip` | 15,269,897 | <https://www.texturecan.com> |

The three `Grass`/`Ground` sets are ambientCG's own export: each ships a
`.tres` whose resource uids carry that site's `acg_` prefix, and a `.mtlx`
naming the maps. `ground_0046` ships its own `ground_0046_description.txt`:

> Name:
> Damaged Concrete with Deep Cracks (Ground 0046)
>
> Description:
> Cracked concrete full of large and small gaps. With the SBSAR file, this
> texture can turn into a dry muddy field or cracky and sandy ground.
>
> Author:
> https://www.texturecan.com

That is what an excavated sandbox's retaining wall should look like, which is
why this set is layer 3 and not one of the others.

**`Ground051` ships no 1K variant.** 2K is the smallest published size, so it is
downscaled from 2048 like the rest are downscaled from 1024; the pack resizes
every map to 512 regardless, so the source size only changes how much detail is
thrown away.

**`NormalGL`, never `NormalDX`.** `triplanar.wgsl`'s whiteout blend adds each
plane's tangent normal to the geometric normal assuming OpenGL's `+Y` up, and
`StandardMaterial::flip_normal_map_y` exists precisely because DirectX maps are
the other one. A DirectX map here inverts the lighting on every slope, which
reads as the terrain being lit from below.

### The sixteen source maps, by SHA-256

| layer | file | bytes | sha256 |
|---|---|---|---|
| 0 | `Grass002_1K-JPG_Color.jpg` | 1,698,232 | `9d6d920f27a38376747e0ac9458ee6449a94d382d2626b3d215d769231a023fd` |
| 0 | `Grass002_1K-JPG_Roughness.jpg` | 844,685 | `3fd5990ed6677e8dde54d4e56e99697a770ccafcedaea74bd538eeb6cef4cf6b` |
| 0 | `Grass002_1K-JPG_NormalGL.jpg` | 2,332,954 | `23ab371ce30f33f54fca8d4c9a0a35d7e70f54ad910e6e7eb5f747b989ab9889` |
| 0 | `Grass002_1K-JPG_AmbientOcclusion.jpg` | 899,701 | `3c15d3e41a49cdc4bc8d6884e9a41394e41ae257fd4aceca9e596adcc5c042c6` |
| 1 | `Ground023_1K-JPG_Color.jpg` | 1,820,307 | `458156780810eefc52a945b55f11d02571042fcd077134b59bb330335f68571a` |
| 1 | `Ground023_1K-JPG_Roughness.jpg` | 746,375 | `639ab71c514098ba7cb1ebd017f1c7a20b7fd0e38a7ac2c5ff55e7bf2ae2f3d7` |
| 1 | `Ground023_1K-JPG_NormalGL.jpg` | 2,521,545 | `7dcdf120d837ccbd48de1a46619fa2b21d2988adf48ce9712a0028448719dc1c` |
| 1 | `Ground023_1K-JPG_AmbientOcclusion.jpg` | 732,506 | `9acacb7c64099c3e570be8a58e331b1cf76a1afb80b4fdf9a0e3c1ad46d05cf3` |
| 2 | `Ground051_2K-JPG_Color.jpg` | 4,069,687 | `49990335422de1938bc504337372100fc7ead50d2370cfde377a1ede634972e5` |
| 2 | `Ground051_2K-JPG_Roughness.jpg` | 1,727,779 | `eb674faff9111ec882dee0e53881f3aa86b1464f8dde92c13a1eac8f5b1364fc` |
| 2 | `Ground051_2K-JPG_NormalGL.jpg` | 8,858,730 | `65466c25b63c7518663e96674f23e8872e5b52975a6615f64b46c248a59ffdd9` |
| 2 | `Ground051_2K-JPG_AmbientOcclusion.jpg` | 2,014,584 | `88d0c9d68bc45f9eb2b89b18ec2abcad7ebf66cb08e3acf254b4e4c78c15aa52` |
| 3 | `ground_0046_color_1k.jpg` | 131,701 | `2b87bff977bcb0413a5fb5e8d383f5b97e3445abf6b39cde1be6898441afa0e0` |
| 3 | `ground_0046_roughness_1k.jpg` | 81,903 | `8af7abf7518f3836657813442df5f04816e9b3c927a324787b428fff0605eacc` |
| 3 | `ground_0046_normal_opengl_1k.png` | 6,494,250 | `961fc6e376324ae1c86e75e75abd5c1270a36d81ca56b8307a77449f6d8eedce` |
| 3 | `ground_0046_ao_1k.jpg` | 236,681 | `63dadc0b500f36ab63fcd7b8f54ce43127c5856acbe28d6c7bc85a71ebaf1ffb` |

`ground_0046_normal_opengl_1k.png` is 6 MB because it is 1024×1024 **16-bit**
RGBA (PNG IHDR bit depth 16, colour type 6). The pack truncates it to 8 bits.

## The pack

| file | RGB | A | layers | `is_srgb` | wgpu format |
|---|---|---|---|---|---|
| `terrain_albedo_roughness_array.png` | colour | roughness | 4 × 512² | `true` | `Rgba8UnormSrgb` |
| `terrain_normal_ao_array.png` | normal (OpenGL) | ambient occlusion | 4 × 512² | `false` | `Rgba8Unorm` |

`Rgba8UnormSrgb` applies the sRGB transfer function to RGB and leaves A linear,
which is exactly right: colour is sRGB-encoded and roughness is linear data.

**512² per layer, not 1024².** Four layers at 1024² is 16 MB of decoded RGBA per
image before mipmaps, decoded on the main thread at startup by
`Image::from_buffer`; at 512² it is 4 MB. The terrain is tiled at 1.5 world
units per tile, so a 512 tile still puts 341 texels on a world unit — more than
a 993×558 canvas can resolve at any distance the player stands.

ImageMagick 7.1.2-29, run from the repository root:

```bash
PBR=/mnt/codex_fs/game_assets/textures/pbr
T=/tmp/terrain_array
OUT=bevy_isomesh/examples/textures
rm -rf "$T"; mkdir -p "$T/src" "$OUT"

unzip -j -o "$PBR/Grass002_1K-JPG.zip" \
  'Grass002_1K-JPG_Color.jpg' 'Grass002_1K-JPG_Roughness.jpg' \
  'Grass002_1K-JPG_NormalGL.jpg' 'Grass002_1K-JPG_AmbientOcclusion.jpg' -d "$T/src"
unzip -j -o "$PBR/Ground023_1K-JPG.zip" \
  'Ground023_1K-JPG_Color.jpg' 'Ground023_1K-JPG_Roughness.jpg' \
  'Ground023_1K-JPG_NormalGL.jpg' 'Ground023_1K-JPG_AmbientOcclusion.jpg' -d "$T/src"
unzip -j -o "$PBR/Ground051_2K-JPG.zip" \
  'Ground051_2K-JPG_Color.jpg' 'Ground051_2K-JPG_Roughness.jpg' \
  'Ground051_2K-JPG_NormalGL.jpg' 'Ground051_2K-JPG_AmbientOcclusion.jpg' -d "$T/src"
unzip -j -o "$PBR/ground_0046_1k_38H6pY.zip" \
  'ground_0046_color_1k.jpg' 'ground_0046_roughness_1k.jpg' \
  'ground_0046_normal_opengl_1k.png' 'ground_0046_ao_1k.jpg' -d "$T/src"

# `-set colorspace sRGB` relabels without transforming; a bare `-colorspace`
# would transform, which would gamma-shift the roughness bytes. It is on the
# colour pack only: the normal pack is direction data loaded with
# `is_srgb: false`. `+profile '*'` drops the ICC profile on the way out --
# nothing downstream reads it.
pack() {
  L=$1; COLOR=$2; ROUGH=$3; NORMAL=$4; AO=$5
  magick "$T/src/$COLOR" -resize '512x512!' -depth 8 -set colorspace sRGB \
    \( "$T/src/$ROUGH" -resize '512x512!' -depth 8 -colorspace Gray \) \
    -alpha off -compose CopyOpacity -composite "$T/ar_$L.png"
  magick "$T/src/$NORMAL" -resize '512x512!' -depth 8 \
    \( "$T/src/$AO" -resize '512x512!' -depth 8 -colorspace Gray \) \
    -alpha off -compose CopyOpacity -composite "$T/na_$L.png"
}

pack 0 Grass002_1K-JPG_Color.jpg  Grass002_1K-JPG_Roughness.jpg \
       Grass002_1K-JPG_NormalGL.jpg  Grass002_1K-JPG_AmbientOcclusion.jpg
pack 1 Ground023_1K-JPG_Color.jpg Ground023_1K-JPG_Roughness.jpg \
       Ground023_1K-JPG_NormalGL.jpg Ground023_1K-JPG_AmbientOcclusion.jpg
pack 2 Ground051_2K-JPG_Color.jpg Ground051_2K-JPG_Roughness.jpg \
       Ground051_2K-JPG_NormalGL.jpg Ground051_2K-JPG_AmbientOcclusion.jpg
pack 3 ground_0046_color_1k.jpg   ground_0046_roughness_1k.jpg \
       ground_0046_normal_opengl_1k.png ground_0046_ao_1k.jpg

# Layer 0 first, because `reinterpret_stacked_2d_as_array` slices top-down.
magick "$T"/ar_0.png "$T"/ar_1.png "$T"/ar_2.png "$T"/ar_3.png -append \
  -define png:color-type=6 -define png:bit-depth=8 +profile '*' \
  "$OUT/terrain_albedo_roughness_array.png"
magick "$T"/na_0.png "$T"/na_1.png "$T"/na_2.png "$T"/na_3.png -append \
  -define png:color-type=6 -define png:bit-depth=8 +profile '*' \
  "$OUT/terrain_normal_ao_array.png"
```

## What the outputs are, measured

```
terrain_albedo_roughness_array.png  2,356,431 bytes  512x2048 depth=8 colourtype=6
terrain_normal_ao_array.png         3,019,521 bytes  512x2048 depth=8 colourtype=6
                                    5,375,952 bytes combined
```

`512x2048` is the gate that matters most and the one a renderer cannot report:
`reinterpret_stacked_2d_as_array(4)` requires `height % 4 == 0` and divides,
so a stack of three or five layers, or a stack built bottom-up, is a silent
mis-index rather than an error. `game_dig`'s
`the_terrain_array_is_four_square_layers_stacked_top_down` asserts
`width == 512`, `height == 4 * width` and
`depth_or_array_layers == 4` after the reinterpretation, on the committed bytes.

The second gate is that the alpha channel really carries the packed map instead
of 255 everywhere — the failure mode of a mis-ordered `-compose CopyOpacity`.
**Per layer**, because a single flat layer averages away against three good ones:

| layer | albedo/roughness A (roughness) | normal/AO A (ambient occlusion) |
|---|---|---|
| 0 grass | 0.363332  σ 0.0557858 | 0.812495  σ 0.0615883 |
| 1 surface dirt | 0.471815  σ 0.0716365 | 0.897955  σ 0.110288 |
| 2 deep dirt | 0.680892  σ 0.0324008 | 0.948841  σ 0.0402447 |
| 3 concrete | 0.738711  σ 0.0460187 | 0.885381  σ 0.104887 |
| whole image | 0.563688  σ 0.161527 | 0.886168  σ 0.097561 |

A standard deviation of `0` on any row would mean that layer's alpha is flat and
its packed map did not make it in. Layer 3's numbers are within 2e-6 of the
`0.738713 / 0.885382` the previous single-set pack measured at 1024², which is
the downscale doing nothing but resampling — the same bytes, the same map, in a
new slot.

The grass roughness mean of 0.36 against concrete's 0.74 is the material
difference showing up as data: blades of grass have specular highlights and a
cracked concrete slab does not.
