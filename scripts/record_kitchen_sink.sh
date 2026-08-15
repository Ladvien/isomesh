#!/usr/bin/env bash
#
# Rebuild `docs/gifs/kitchen-sink.gif` — the README's banner, several examples
# running at once.
#
#   scripts/record_kitchen_sink.sh [output.gif]
#
# The banner existed and the recipe did not, exactly as for the individual GIFs
# (E-215). This is that recipe: it records each panel through the same capture
# rig `record_gif.sh` uses, then tiles the frame sequences with ffmpeg's `xstack`.
#
# # Why the panels are recorded rather than cropped from one run
#
# Each panel is a different example and they cannot run in one process. So each
# is captured separately at panel resolution and the sequences are stacked
# frame-for-frame, which also means every panel has the same frame count and the
# loop closes cleanly.
#
# # Panel choice
#
# Eight, in a 4x2 grid, chosen to span what the crate does rather than to look
# busy: authoring a field, meshing it, the hard topology case, the two accuracy
# comparisons, and three game-shaped uses. `sdf_authoring` is first among the
# meshing panels on purpose — E-216's point is that nothing showed the SDF as a
# medium rather than an input.

set -euo pipefail
cd "$(dirname "$0")/.."

OUTPUT="${1:-docs/gifs/kitchen-sink.gif}"

: "${DISPLAY:=:0}"
: "${PANEL_W:=320}"
: "${PANEL_H:=180}"
: "${FPS:=20}"
export DISPLAY

# Each entry: example name, then any example-specific environment.
PANELS=(
    "sdf_authoring|ISOMESH_SPIN=0.010"
    "game_showcase|ISOMESH_SPIN=0.004"
    "subgrid_features|ISOMESH_SPIN=0.012"
    "marching_cubes_tunnel|ISOMESH_SPIN=0.012"
    "surface_nets_vs_marching_cubes|ISOMESH_SPIN=0.012"
    "game_dig|"
    "game_walk|"
    "game_destruction|"
)

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Capture at 16:9 and downscale, rather than asking for a 320x180 window: a Bevy
# window that small puts the HUD over the whole frame and the geometry vanishes.
export ISOMESH_WINDOW="${ISOMESH_WINDOW:-1280x720}"
export ISOMESH_CAPTURE_FRAMES="${ISOMESH_CAPTURE_FRAMES:-80}"
export ISOMESH_CAPTURE_EVERY="${ISOMESH_CAPTURE_EVERY:-2}"
export ISOMESH_CAPTURE_SETTLE="${ISOMESH_CAPTURE_SETTLE:-60}"

index=0
for entry in "${PANELS[@]}"; do
    example="${entry%%|*}"
    extra="${entry#*|}"
    dir="$WORK/panel$index"
    mkdir -p "$dir"
    echo "-- [$((index + 1))/${#PANELS[@]}] $example"
    (
        export ISOMESH_CAPTURE="$dir"
        # shellcheck disable=SC2086 -- extra is a deliberate word list of VAR=VAL
        cd bevy_isomesh && env $extra cargo run --example "$example" --release >/dev/null 2>&1
    ) || true
    count="$(find "$dir" -name 'frame_*.png' | wc -l)"
    if [ "$count" -eq 0 ]; then
        echo "::error::$example captured no frames" >&2
        exit 1
    fi
    index=$((index + 1))
done

# Shortest panel decides the length, so every panel is present in every frame.
SHORTEST=1000000
for i in $(seq 0 $((index - 1))); do
    n="$(find "$WORK/panel$i" -name 'frame_*.png' | wc -l)"
    [ "$n" -lt "$SHORTEST" ] && SHORTEST="$n"
done
echo "-- $index panels, $SHORTEST frames each"

INPUTS=()
SCALES=""
STACK=""
for i in $(seq 0 $((index - 1))); do
    INPUTS+=(-framerate "$FPS" -start_number 0 -i "$WORK/panel$i/frame_%04d.png")
    SCALES="$SCALES[$i:v]scale=$PANEL_W:$PANEL_H:flags=lanczos,setsar=1[p$i];"
    STACK="$STACK[p$i]"
done

LAYOUT="0_0|w0_0|w0+w1_0|w0+w1+w2_0|0_h0|w0_h0|w0+w1_h0|w0+w1+w2_h0"

ffmpeg -hide_banner -loglevel error -y "${INPUTS[@]}" \
    -filter_complex "${SCALES}${STACK}xstack=inputs=$index:layout=$LAYOUT[tiled];[tiled]palettegen=max_colors=128:stats_mode=diff[pal]" \
    -map "[pal]" -frames:v 1 "$WORK/palette.png"

ffmpeg -hide_banner -loglevel error -y "${INPUTS[@]}" -i "$WORK/palette.png" \
    -filter_complex "${SCALES}${STACK}xstack=inputs=$index:layout=$LAYOUT[tiled];[tiled][$index:v]paletteuse=dither=bayer:bayer_scale=3[out]" \
    -map "[out]" -frames:v "$SHORTEST" "$OUTPUT"

SIZE="$(stat -c%s "$OUTPUT" 2>/dev/null || stat -f%z "$OUTPUT")"
awk "BEGIN{printf \"-- wrote $OUTPUT (%.2f MB)\\n\", $SIZE/1048576}"
