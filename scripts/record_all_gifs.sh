#!/usr/bin/env bash
#
# Record every committed GIF, from one command.
#
# `record_gif.sh` records **one** clip and takes its parameters from the
# environment. Which parameters each example wants is not obvious and is not
# written down anywhere else: `game_dig` is a static landscape without
# `ISOMESH_AUTOCARVE`, `game_destruction` fires 24 shots over 21 seconds and the
# default capture window catches three, and `marching_cubes_interior` loops once
# every 8.3 seconds so a short window photographs an arbitrary arc of it.
#
# That knowledge lived in a subagent's inventory and one earlier session's shell
# history. This file is where it lives now.
#
#   scripts/record_all_gifs.sh              every clip
#   scripts/record_all_gifs.sh dual_ qef_   only clips whose name matches
#
# Needs an X display. The capture is a GPU readback through Bevy's screenshot
# path rather than a screen scrape, so the window may be unmapped and no window
# manager is required -- but the window must exist, so `DISPLAY` must resolve.
# On this machine that is `:1`, the Xwayland server inside the Hyprland session.

set -uo pipefail

cd "$(dirname "$0")/.."

: "${DISPLAY:=:1}"
export DISPLAY

OUT=docs/gifs
mkdir -p "$OUT"

# Every clip: output stem, example, then the environment it needs.
#
# `nohud` on anything meant to be *looked at* rather than read -- the harness
# added that flag for exactly this. It stays **off** where the HUD is the
# content, which is `game_lod_flyover` (seam counts per side) and
# `resolution_plot` (the fit).
CLIPS=(
    # -- the seven whose sweep code was written for this and never used --------
    "dual-contouring-vs-surface-nets|dual_contouring_cube|ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=60 ISOMESH_CAPTURE_EVERY=2"
    "sharp-features-lambda-sweep|sharp_features|ISOMESH_VIEW=nohud ISOMESH_FIELD=3 ISOMESH_CAPTURE_FRAMES=96 ISOMESH_CAPTURE_EVERY=2"
    "qef-clamp-self-intersections|qef_clamp|ISOMESH_FIELD=2 ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "precision-f32-tears|precision_f32_vs_f64|ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "manifold-check-resolution|manifold_check|ISOMESH_ALGORITHM=sn ISOMESH_FIELD=5 ISOMESH_CAPTURE_FRAMES=72 ISOMESH_CAPTURE_EVERY=2"
    "ambiguous-faces-are-rare|marching_cubes_ambiguity|ISOMESH_FIELD=0 ISOMESH_CAPTURE_FRAMES=60 ISOMESH_CAPTURE_EVERY=2"
    "undo-is-a-refold|game_editor|ISOMESH_CAPTURE_FRAMES=64 ISOMESH_CAPTURE_EVERY=2"

    # -- self-animating, no flags needed --------------------------------------
    # One full sweep is ~8.3 s; 80x6 ticks is ~8 s at 60 Hz.
    "the-interior-decider-sweep|marching_cubes_interior|ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=6"
    "a-boolean-remeshed-every-frame|game_csg_props|ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=4"
    # Fires every 0.9 s up to 24 shots. 100x8 ticks is ~13 s, about 14 shots.
    "the-debris-is-the-boolean|game_destruction|ISOMESH_TARGET=shell ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=100 ISOMESH_CAPTURE_EVERY=8"
    "lod-flyover|game_lod_flyover|ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=3"
    "gpu-resident-mesh-shader|gpu_mesh_shader|ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=3"
    # The 89^3 points dominate, so it fills in slowly. Long stride, few frames.
    "the-fit-drawing-itself|resolution_plot|ISOMESH_CAPTURE_FRAMES=60 ISOMESH_CAPTURE_EVERY=4"

    # -- need a flag to be worth filming --------------------------------------
    "digging-a-tunnel|game_dig|ISOMESH_AUTOCARVE=60 ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "paint-that-survives-the-wall|game_paint|ISOMESH_AUTOPAINT=60 ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"

    # -- the ten that already existed, re-recorded at the current commit ------
    "flying-through-the-rock|game_showcase|ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=3"
    "walking-the-seams|game_walk|ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=3"
    "terrain-streaming|game_terrain_stream|ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=3"
    "building-a-field|sdf_authoring|ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "the-tunnel-meshed-as-a-tunnel|marching_cubes_tunnel|ISOMESH_SPIN=0.012 ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "subgrid-letters-thinner-than-a-voxel|subgrid_features|ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "surface-nets-vs-marching-cubes-box|surface_nets_vs_marching_cubes|ISOMESH_FIELD=2 ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "surface-nets-vs-marching-cubes-gyroid|surface_nets_vs_marching_cubes|ISOMESH_FIELD=3 ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "marching-cubes-sphere-resolution-sweep|marching_cubes_sphere|ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
)

# The filters are this script's own arguments, captured once. `want` used to
# take them as trailing parameters and test `$#` for emptiness, which is always
# false inside a function that was also handed the stem -- so "no filter" meant
# "match nothing" and a full run recorded zero clips.
FILTERS=("$@")

want() {
    [ "${#FILTERS[@]}" -eq 0 ] && return 0
    local stem=$1 pattern
    for pattern in "${FILTERS[@]}"; do
        case "$stem" in *"$pattern"*) return 0 ;; esac
    done
    return 1
}

ok=0
skipped=0
failed=()
for clip in "${CLIPS[@]}"; do
    IFS='|' read -r stem example environment <<<"$clip"
    if ! want "$stem"; then
        skipped=$((skipped + 1))
        continue
    fi
    printf '\n\033[1m== %s  (%s)\033[0m\n' "$stem" "$example"
    # shellcheck disable=SC2086
    if env $environment ./scripts/record_gif.sh "$example" "$OUT/$stem.gif"; then
        ok=$((ok + 1))
    else
        failed+=("$stem")
    fi
done

printf '\n%s recorded' "$ok"
[ "$skipped" -gt 0 ] && printf ', %s skipped' "$skipped"
if [ "${#failed[@]}" -gt 0 ]; then
    printf ', \033[31m%s failed\033[0m: %s\n' "${#failed[@]}" "${failed[*]}"
    exit 1
fi
printf '\n'
# The kitchen-sink banner is eight of these panels stacked, so it goes last and
# has its own script.
printf 'now rebuild the banner: scripts/record_kitchen_sink.sh\n'
