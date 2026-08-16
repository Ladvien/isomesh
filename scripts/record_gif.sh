#!/usr/bin/env bash
#
# Record a Bevy example to a GIF.
#
#   scripts/record_gif.sh <example> <output.gif> [ISOMESH_* overrides as env]
#
# Nine GIFs existed before this script and the method that made them was written
# down nowhere -- grepping `docs/` and `scripts/` for `ffmpeg`, `x11grab` and
# `palettegen` returned nothing. This is E-215: the tenth is a command.
#
# # It does not scrape the screen
#
# The obvious approach is `ffmpeg -f x11grab`, and it is the wrong one. The
# examples already carry a capture rig -- `ISOMESH_CAPTURE` writes a numbered
# frame per N ticks through Bevy's own screenshot path, which reads back from the
# GPU. That works over a window the compositor never mapped, needs no window
# manager, and cannot catch another window passing over the top. This script
# drives that rig and assembles its output.
#
# # Two things that are measured rather than chosen
#
# `ISOMESH_WINDOW` sets the size, and it only works because `size_window` runs in
# `Update` and re-applies across the frames in which the window is created -- at
# `PreStartup` the OS window does not exist yet and the request is silently lost
# (E-214, FINDINGS M-235). If a capture comes back the wrong shape, that is the
# first thing to check.
#
# The GIF is built in two passes. A single pass quantises against a fixed palette
# and bands badly on shaded 3D; `palettegen`/`paletteuse` builds the palette from
# the footage. The committed GIFs run 0.7-4.8 MB and this warns outside that.

set -euo pipefail
cd "$(dirname "$0")/.."

if [ "$#" -lt 2 ]; then
    echo "usage: $0 <example> <output.gif>" >&2
    echo "  env: ISOMESH_WINDOW=1280x720 ISOMESH_SPIN=0.012 ISOMESH_VIEW=nohud" >&2
    echo "       ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2 ISOMESH_CAPTURE_SETTLE=60" >&2
    echo "       FPS=20 WIDTH=900 DISPLAY=:1" >&2
    exit 2
fi

EXAMPLE="$1"
OUTPUT="$2"

: "${DISPLAY:=:0}"
: "${ISOMESH_WINDOW:=1280x720}"
: "${ISOMESH_CAPTURE_FRAMES:=80}"
: "${ISOMESH_CAPTURE_EVERY:=2}"
# Comfortably past SIZE_WINDOW_FRAMES, so frame zero is taken after the window
# has stopped changing size.
: "${ISOMESH_CAPTURE_SETTLE:=60}"
: "${FPS:=20}"
: "${WIDTH:=900}"
# Palette size and dither.
#
# These used to be hard-coded at `128` and `bayer:bayer_scale=3`, and that was
# the wrong trade for almost every clip here. `bayer` is an *ordered* dither:
# cheap, very compressible, and it lays a fixed crosshatch over the whole frame
# that reads as chunky on smooth shading -- which is most of what this project
# renders. Side by side on the terrain flythrough the difference is not subtle.
#
# `sierra2_4a` diffuses the error instead. It looks far better and costs
# roughly 2-3x the bytes, because error diffusion decorrelates neighbouring
# frames and defeats the inter-frame compression a GIF relies on. So the good
# setting is the default, and the handful of clips where full-frame motion makes
# that unaffordable ask for the cheap one by name.
: "${COLORS:=256}"
: "${DITHER:=sierra2_4a}"
export DISPLAY ISOMESH_WINDOW ISOMESH_CAPTURE_FRAMES ISOMESH_CAPTURE_EVERY ISOMESH_CAPTURE_SETTLE

for tool in ffmpeg cargo; do
    command -v "$tool" >/dev/null || {
        echo "::error::$tool is not installed" >&2
        exit 1
    }
done

FRAMES="$(mktemp -d)"
trap 'rm -rf "$FRAMES"' EXIT
export ISOMESH_CAPTURE="$FRAMES"

echo "-- recording $EXAMPLE at $ISOMESH_WINDOW on DISPLAY=$DISPLAY"
# Release, always. A debug build meshes 37-62x slower (M-152) and the capture
# would photograph a mid-extraction frame.
( cd bevy_isomesh && cargo run --example "$EXAMPLE" --release >/dev/null 2>&1 ) || true

COUNT="$(find "$FRAMES" -name 'frame_*.png' | wc -l)"
if [ "$COUNT" -eq 0 ]; then
    echo "::error::no frames captured -- is DISPLAY=$DISPLAY reachable?" >&2
    exit 1
fi
echo "-- $COUNT frames"

PALETTE="$FRAMES/palette.png"
ffmpeg -hide_banner -loglevel error -y -framerate "$FPS" -i "$FRAMES/frame_%04d.png" \
    -vf "scale=$WIDTH:-1:flags=lanczos,palettegen=max_colors=$COLORS:stats_mode=diff" "$PALETTE"
ffmpeg -hide_banner -loglevel error -y -framerate "$FPS" -i "$FRAMES/frame_%04d.png" -i "$PALETTE" \
    -lavfi "scale=$WIDTH:-1:flags=lanczos[s];[s][1:v]paletteuse=dither=$DITHER" \
    "$OUTPUT"

SIZE="$(stat -c%s "$OUTPUT" 2>/dev/null || stat -f%z "$OUTPUT")"
MB="$(awk "BEGIN{printf \"%.2f\", $SIZE/1048576}")"
echo "-- wrote $OUTPUT ($MB MB)"
# The committed GIFs sit in 0.7-4.8 MB. Outside that is not an error -- a shorter
# clip is legitimately smaller -- but it is worth seeing before committing.
awk "BEGIN{exit !($SIZE > 5033164)}" && echo "::warning::$MB MB is above the 4.8 MB the committed GIFs sit within"
awk "BEGIN{exit !($SIZE < 734003)}" && echo "::warning::$MB MB is below the 0.7 MB the committed GIFs sit within"
exit 0
