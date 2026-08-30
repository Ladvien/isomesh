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
# **`ISOMESH_FIELD` is a zero-based index into each example's own list, and the
# lists differ.** `surface_nets_vs_marching_cubes` is [sphere, box_exact, torus,
# csg_difference, gyroid]; `sharp_features` is [box_exact, csg_difference,
# gyroid, fbm_terrain]; `qef_clamp` is [gyroid, torus, fbm_terrain,
# csg_difference]; `manifold_check` uses the crate's reference-field order.
# Guessing cost three clips filed under names their contents did not match.
#
# # `nohud` only where nothing on screen needs naming
#
# The harness added that flag for "a GIF meant to be looked at rather than
# read", and a first pass applied it to every clip that was not obviously a
# numbers demo. That was wrong for a whole class: **on a side-by-side comparison
# the label is the content.** `dual_contouring_cube` puts Surface Nets and Dual
# Contouring next to each other in two shades, and without the HUD a viewer
# cannot tell which is which -- so the clip shows two objects and proves nothing.
# Same for both `surface_nets_vs_marching_cubes` clips, for `subgrid_features`
# (the letter thickness is the sweep) and `sharp_features` (lambda is).
#
# So `nohud` is for the *worlds* -- showcase, terrain, dig, paint, destruction,
# csg props, mesh shader -- where the picture describes itself. Anything that
# compares, sweeps or counts keeps its HUD.
#
# # Two calibrations that a first pass got wrong, both measured
#
# **A scripted sequence must outlast the settle.** `record_gif.sh` waits 60 ticks
# before frame zero, and `ISOMESH_AUTOCARVE=60` carves one per frame -- so the
# whole tunnel was dug before the recording began and the clip was a still of the
# end state. `game_paint` and `resolution_plot` failed the same way. The window is
# `SETTLE + FRAMES x EVERY` ticks, and the animation has to cover it.
#
# **A wide stride is what makes a GIF enormous.** `ISOMESH_CAPTURE_EVERY=3` on the
# flying demos put so much motion between frames that the palette could not
# compress them: `terrain-streaming` came back at **24 MB** against a 2.7 MB
# predecessor. Stride 2 and a narrower `WIDTH` fix it. More frames of a smaller
# picture beat fewer frames of a bigger one, every time.
CLIPS=(
    # -- the seven whose sweep code was written for this and never used --------
    "dual-contouring-vs-surface-nets|dual_contouring_cube|ISOMESH_CAPTURE_FRAMES=60 ISOMESH_CAPTURE_EVERY=2"
    # `csg_difference`, after trying the other two. `box_exact` is grid-aligned
    # enough that the sweep barely moves a pixel, and `gyroid` fills the frame at
    # the harness's default camera radius of 7 -- you end up inside the surface.
    # A box with a sphere bitten out of it has both a convex rim and a concave
    # one, at a scale the default camera frames, and both visibly round over.
    "sharp-features-lambda-sweep|sharp_features|ISOMESH_FIELD=1 ISOMESH_CAPTURE_FRAMES=48 ISOMESH_CAPTURE_EVERY=2 WIDTH=700 FPS=8"
    "qef-clamp-self-intersections|qef_clamp|ISOMESH_FIELD=0 ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
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
    "lod-flyover|game_lod_flyover|ISOMESH_CAPTURE_FRAMES=64 ISOMESH_CAPTURE_EVERY=2 WIDTH=800"
    "gpu-resident-mesh-shader|gpu_mesh_shader|ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=3"

    # -- need a flag to be worth filming --------------------------------------
    "digging-a-tunnel|game_dig|ISOMESH_AUTOCARVE=16 ISOMESH_AUTOCARVE_EVERY=3 ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=60 ISOMESH_CAPTURE_EVERY=2 FPS=10"
    "paint-that-survives-the-wall|game_paint|ISOMESH_AUTOPAINT=38 ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=60 ISOMESH_CAPTURE_EVERY=2 FPS=10"

    # -- the ten that already existed, re-recorded at the current commit ------
    # Three knobs, and each fixes something a first pass got wrong. **Speed
    # 3.0** rather than the interactive default of 5: at 5 a chunk goes from
    # appearing to filling the frame in about a second and the clip reads as
    # flashing. **A 48-metre stream radius** rather than 34, so arrivals happen
    # in the far distance where they belong rather than mid-frame. And a **220
    # tick settle**, because slowing the flight also means it has covered less
    # ground by the time recording starts -- at the default settle the camera is
    # still over open plain and half the frame is sky.
    "flying-through-the-rock|game_showcase|ISOMESH_VIEW=nohud ISOMESH_SPEED=3.0 ISOMESH_STREAM_VIEW=48 ISOMESH_CAPTURE_SETTLE=220 ISOMESH_CAPTURE_FRAMES=48 ISOMESH_CAPTURE_EVERY=2 WIDTH=1000"
    "walking-the-seams|game_walk|ISOMESH_CAPTURE_FRAMES=48 ISOMESH_CAPTURE_EVERY=2 WIDTH=640"
    # The narrowest clip here, and it has to be: endless noisy terrain flying
    # past is the worst case a GIF palette can be handed. Every pixel changes
    # every frame and none of it repeats, so this is the one demo where the
    # picture has to shrink rather than the clip get shorter.
    "terrain-streaming|game_terrain_stream|ISOMESH_CAPTURE_FRAMES=40 ISOMESH_CAPTURE_EVERY=2 WIDTH=540"
    "building-a-field|sdf_authoring|ISOMESH_CAPTURE_FRAMES=72 ISOMESH_CAPTURE_EVERY=2"
    "the-tunnel-meshed-as-a-tunnel|marching_cubes_tunnel|ISOMESH_SPIN=0.012 ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "subgrid-letters-thinner-than-a-voxel|subgrid_features|ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "surface-nets-vs-marching-cubes-box|surface_nets_vs_marching_cubes|ISOMESH_FIELD=1 ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "surface-nets-vs-marching-cubes-gyroid|surface_nets_vs_marching_cubes|ISOMESH_FIELD=4 ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "marching-cubes-sphere-resolution-sweep|marching_cubes_sphere|ISOMESH_VIEW=nohud ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"

    # -- Phase 21 -------------------------------------------------------------
    # `game_edit_tape_trim` runs the whole 1,571-re-mesh ablation before the
    # window opens, so the settle is spent on a splash rather than on the scene;
    # 150 frames at stride 2 covers the fifteen-segment tour, and the tour CUTS
    # between chunks rather than panning, which is why 150 frames costs the same
    # 0.92 MB that 120 did.
    "the-tape-you-keep-is-twenty-times-too-big|game_edit_tape_trim|ISOMESH_CAPTURE_FRAMES=150 ISOMESH_CAPTURE_EVERY=2 FPS=10"
    "mirrored-is-not-the-same-mesh|game_mirror_dedup|ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
    "where-the-root-falls-decides-the-gain|shifted_linear_root|ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"

    # -- Phase 27 -------------------------------------------------------------
    # All five are self-driving under capture, so every row is the demo's own
    # walk. The HUD stays on for all five because on a census or a sweep the
    # label IS the content -- hyperdeterminant compares a live panel to
    # P-130's partition on the same screen, anisotropic_metric's two panels are
    # compared by the numbers beside them. hyperdeterminant's 96 frames cover
    # its eight-field walk at twelve a field; tpms's 90 cover the nine-state
    # cycle at ten a state; cave's 82 cover the isovalue ladder across the
    # transition; intrinsic's 80 give four full before/after alternations at
    # ten frames an arm.
    "the-hyperdeterminant-in-every-cell|hyperdeterminant_cells|ISOMESH_CAPTURE_FRAMES=96 ISOMESH_CAPTURE_EVERY=2"
    "three-periodic-surfaces-with-a-known-topology|tpms_euler|ISOMESH_CAPTURE_FRAMES=90 ISOMESH_CAPTURE_EVERY=2"
    "the-metric-that-costs-more-than-it-saves|anisotropic_metric|ISOMESH_CAPTURE_FRAMES=96 ISOMESH_CAPTURE_EVERY=2"
    "where-the-caves-join-up|cave_percolation|ISOMESH_FIELD=1 ISOMESH_CAPTURE_FRAMES=82 ISOMESH_CAPTURE_EVERY=2"
    "fifteen-thousand-flips-that-move-nothing|intrinsic_flips|ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2"
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
