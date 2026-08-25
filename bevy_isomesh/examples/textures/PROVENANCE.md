# `ground_0046` — where these two PNGs came from

`game_dig` textures its rock with a triplanar blend of the two files in this
directory. They are not the artist's files: they are four of that set's maps
packed into two RGBA images, because WebGL2 counts bound textures and
`StandardMaterial` already spends six of them. This file is the record of what
was packed and how, so the pack can be redone or a different set swapped in
without guessing.

The bytes are committed rather than fetched at build time. The Open SciVis
volumes `.gitignore` excludes (M-006) have a public URL to fetch from; this set
does not, and the repository already tracks comparable binaries — the largest
tracked file is `docs/gifs/flying-through-the-rock.gif` at 7.4 MB.

## Source

`/mnt/codex_fs/game_assets/textures/pbr/ground_0046_1k_38H6pY.zip`
(15,269,897 bytes), from <https://www.texturecan.com>.

The set's own `ground_0046_description.txt`:

> Name:
> Damaged Concrete with Deep Cracks (Ground 0046)
>
> Description:
> Cracked concrete full of large and small gaps. With the SBSAR file, this
> texture can turn into a dry muddy field or cracky and sandy ground.
>
> Author:
> https://www.texturecan.com

Its `ground_0046_keywords.txt`:

> ground, concrete, cracks, cracked, dry, drought, concrete, gaps, damaged,
> grey, gray

That is what an excavated sandbox should look like, which is why this set and
not one of the grass or gravel ones.

`normal_opengl`, not `normal_directx`: Bevy's normal-map convention is `+Y` up,
and `StandardMaterial::flip_normal_map_y` exists precisely because DirectX maps
are the other one.

### The four source maps, by SHA-256

| file | bytes | sha256 |
|---|---|---|
| `ground_0046_color_1k.jpg` | 131,701 | `2b87bff977bcb0413a5fb5e8d383f5b97e3445abf6b39cde1be6898441afa0e0` |
| `ground_0046_roughness_1k.jpg` | 81,903 | `8af7abf7518f3836657813442df5f04816e9b3c927a324787b428fff0605eacc` |
| `ground_0046_ao_1k.jpg` | 236,681 | `63dadc0b500f36ab63fcd7b8f54ce43127c5856acbe28d6c7bc85a71ebaf1ffb` |
| `ground_0046_normal_opengl_1k.png` | 6,494,250 | `961fc6e376324ae1c86e75e75abd5c1270a36d81ca56b8307a77449f6d8eedce` |

The normal map is 6 MB because it is 1024×1024 **16-bit** RGBA (PNG IHDR bit
depth 16, colour type 6). The pack truncates it to 8 bits, which is why the
output is under half the size.

## The pack

| file | RGB | A | `is_srgb` | wgpu format |
|---|---|---|---|---|
| `ground_0046_albedo_roughness.png` | colour | roughness | `true` | `Rgba8UnormSrgb` |
| `ground_0046_normal_ao.png` | normal (OpenGL) | ambient occlusion | `false` | `Rgba8Unorm` |

`Rgba8UnormSrgb` applies the sRGB transfer function to RGB and leaves A linear,
which is exactly right: colour is sRGB-encoded and roughness is linear data.

ImageMagick 7.1.2-29, run from the repository root:

```bash
mkdir -p /tmp/g46 bevy_isomesh/examples/textures
unzip -j -o /mnt/codex_fs/game_assets/textures/pbr/ground_0046_1k_38H6pY.zip \
    'ground_0046_color_1k.jpg' 'ground_0046_roughness_1k.jpg' \
    'ground_0046_ao_1k.jpg' 'ground_0046_normal_opengl_1k.png' -d /tmp/g46

# `-set colorspace sRGB` relabels without transforming; a bare `-colorspace`
# would transform, which would gamma-shift the roughness bytes. `+profile '*'`
# drops the ICC profile: nothing downstream reads it.
magick /tmp/g46/ground_0046_color_1k.jpg -depth 8 -set colorspace sRGB \
    \( /tmp/g46/ground_0046_roughness_1k.jpg -depth 8 -colorspace Gray \) \
    -alpha off -compose CopyOpacity -composite \
    -define png:color-type=6 -define png:bit-depth=8 +profile '*' \
    bevy_isomesh/examples/textures/ground_0046_albedo_roughness.png

# No `-set colorspace` on this one: the normal map is direction data, it is
# loaded with `is_srgb: false`, and the only wanted change is 16-bit to 8-bit.
magick /tmp/g46/ground_0046_normal_opengl_1k.png -depth 8 \
    \( /tmp/g46/ground_0046_ao_1k.jpg -depth 8 -colorspace Gray \) \
    -alpha off -compose CopyOpacity -composite \
    -define png:color-type=6 -define png:bit-depth=8 +profile '*' \
    bevy_isomesh/examples/textures/ground_0046_normal_ao.png
```

## What the outputs are, measured

```
ground_0046_albedo_roughness.png  1,378,415 bytes  1024x1024 depth=8 colourtype=6
ground_0046_normal_ao.png         2,950,213 bytes  1024x1024 depth=8 colourtype=6
                                  4,328,628 bytes combined
```

The gate that matters is not the size, it is that the alpha channel really
carries the packed map instead of 255 everywhere — the failure mode of a
mis-ordered `-compose CopyOpacity`. Alpha mean and standard deviation:

```
ground_0046_albedo_roughness.png   0.738713  0.0484046   (roughness)
ground_0046_normal_ao.png          0.885382  0.12074     (ambient occlusion)
```

A standard deviation of `0` there would mean the alpha is flat and the packed
map did not make it in.
