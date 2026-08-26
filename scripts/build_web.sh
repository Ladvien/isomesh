#!/usr/bin/env bash
#
# Build the whole published site into web/dist: three WebAssembly demos and the
# rendered prose.
#
#   scripts/build_web.sh
#
# This is the single entry point. `scripts/build_site.py` renders the prose and
# can be run on its own while iterating on it; this script owns the clean and the
# wasm half, and calls that one at the end so both land in the same tree.
#
# # Two things that are not obvious
#
# **The `wasm-bindgen` CLI version has to equal the `wasm-bindgen` crate version
# in the lockfile.** A CLI that differs emits glue for a different ABI and the
# module fails to instantiate in the browser with an error that names neither
# tool. So the version is read out of `bevy_isomesh/Cargo.lock` rather than
# written here, and a mismatch is a hard stop with the exact `cargo install` line
# to fix it.
#
# **The demo list is an array**, so a tenth playable demo is one line here plus
# one entry in `web/play.html`'s allow-list -- and `scripts/doc_facts.sh` checks
# those two against each other, because a module built but not allow-listed is
# 36 MB nothing can reach and one allow-listed but not built is a link to a 404.
# Only these nine are built for the web; the other examples are untouched and
# stay native-only.
#
# **`isomesh_web` is not one of them.** It is the front page's own module -- the
# core crate with a hand-written WebGL2 renderer instead of Bevy, ~130 KB
# instead of 36 MB -- and it is built here because this is the entry point, not
# because it is a Bevy example.
#
# `PYTHON` exists because the renderer needs two pip packages and a PEP 668
# distribution will not install them into the system interpreter. Point it at a
# venv's python: `PYTHON=~/.venvs/isomesh/bin/python scripts/build_web.sh`.

set -euo pipefail
cd "$(dirname "$0")/.."

# `game_dig` first because it is the flagship: it is the demo the front page and
# both prose sites lead with, so a break in it is the break that matters most and
# should surface before twenty minutes of the rest. Then cheap and visual, so a
# break in the wasm build in general still surfaces on a two-minute module rather
# than at the end. The three Phase 21 demos are last because they are the ones
# with a cross-check to run.
DEMOS=(
    game_dig
    quickstart
    marching_cubes_tunnel
    dual_contouring_cube
    surface_nets_vs_marching_cubes
    game_showcase
    game_mirror_dedup
    game_edit_tape_trim
    shifted_linear_root
)
OUT=web/dist
PYTHON="${PYTHON:-python3}"

# There is deliberately no `CNAME` file and no custom domain.
#
# `isomesh.ladvien.com` was the custom domain for a while and never got a TLS
# certificate: `gh api repos/Ladvien/isomesh/pages` reported `https_certificate`
# absent while `pages/health` reported `is_https_eligible: true`, `is_valid:
# true`, `caa_error: null` and `https_error: peer_failed_verification` -- GitHub
# agreeing the certificate should exist and no certificate existing. DNS was
# clean and independently verified: a CNAME to `ladvien.github.io`, unproxied, on
# Route 53, with `github.io`'s CAA permitting `letsencrypt.org`.
#
# The consequence was not cosmetic. **WebGPU is a secure-context-only API**, and
# `https://ladvien.github.io/isomesh/` 301-redirected to
# `http://isomesh.ladvien.com/`, so `navigator.gpu` was `undefined` and
# `web/play.html`'s gate fired in every browser -- Chrome on a desktop included.
# Every demo on this site needs WebGPU, so the pretty URL cost the whole site.
#
# The site now serves at `https://ladvien.github.io/isomesh/` under GitHub's own
# `*.github.io` certificate, which is valid. Do not re-add the domain here: a
# `CNAME` file in the artifact *re-sets* the custom domain on every deploy, so a
# line here silently undoes the removal on the next push.
WANT=$(awk '/^name = "wasm-bindgen"$/{getline; gsub(/[^0-9.]/,""); print; exit}' \
    bevy_isomesh/Cargo.lock)
if [ -z "$WANT" ]; then
    echo "cannot read the wasm-bindgen version from bevy_isomesh/Cargo.lock" >&2
    exit 1
fi

if ! command -v wasm-bindgen >/dev/null 2>&1 ||
    ! wasm-bindgen --version 2>/dev/null | grep -qF "$WANT"; then
    echo "wasm-bindgen $WANT is required and is not what is on PATH." >&2
    echo "  cargo install wasm-bindgen-cli --version $WANT --locked" >&2
    exit 1
fi

# `#[unsafe(no_mangle)]` is the only `unsafe` token `isomesh_web` is allowed, and
# this is what enforces it. The crate carries `unsafe_code = "allow"` because
# edition 2024 has no way to emit a wasm export without that attribute, so the
# lint cannot draw the line and this does: an `unsafe` block, `fn`, `impl` or
# `trait` in that crate stops the build here rather than shipping.
#
# It sits with the other refusals, above `rm -rf "$OUT"`, and deliberately: a
# gate that fires *after* the clean has already destroyed the previous build, so
# a local iteration loop pays 432 MB and eight minutes for a typo. Every check
# that can fail without compiling anything belongs on this side of that line.
if grep -rEn 'unsafe *(\{|fn |impl |trait )' isomesh_web/src; then
    echo "isomesh_web may use #[unsafe(no_mangle)] and nothing else -- see its [lints.rust]" >&2
    exit 1
fi

rustup target add wasm32-unknown-unknown

rm -rf "$OUT"
mkdir -p "$OUT/play/pkg"

echo "==> isomesh_web (the front page's module)"
(cd isomesh_web && cargo build --release --target wasm32-unknown-unknown)
# Copied rather than piped through `wasm-opt`, and that is a measurement rather
# than an omission. At M-362 `-Oz` took this module from 133,115 to 112,642 bytes
# raw and from 49,837 to **50,493** gzipped: 15.4% smaller on disk and 1.3%
# *bigger* on the wire, which is the number a front page actually pays. The Bevy
# loop below was measured the same way and is worse -- 37,409,471 to 29,094,422
# raw, 8,729,859 to **9,314,752** gzipped, 6.7% the wrong way, at 23 s per
# module -- so there is no `wasm-opt` anywhere in this script and no `binaryen`
# in CI. Each crate's own `[profile.release]` is doing this work already.
cp isomesh_web/target/wasm32-unknown-unknown/release/isomesh_web.wasm \
    "$OUT/isomesh_web.wasm"

for demo in "${DEMOS[@]}"; do
    echo "==> $demo"
    (cd bevy_isomesh &&
        cargo build --profile wasm-release --target wasm32-unknown-unknown \
            --example "$demo")
    wasm-bindgen --target web --no-typescript \
        --remove-name-section --remove-producers-section \
        --out-dir "$OUT/play/pkg/$demo" --out-name "$demo" \
        "bevy_isomesh/target/wasm32-unknown-unknown/wasm-release/examples/$demo.wasm"
done

echo "==> prose and assets"
"$PYTHON" scripts/build_site.py

echo
echo "wasm modules:"
for demo in "${DEMOS[@]}"; do
    du -h "$OUT/play/pkg/$demo/${demo}_bg.wasm" | sed 's/^/  /'
done
du -h "$OUT/isomesh_web.wasm" | sed 's/^/  /'
echo "site total:"
du -sh "$OUT" | sed 's/^/  /'
