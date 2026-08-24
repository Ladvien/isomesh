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
# **The demo list is an array**, so a fourth playable demo is one line here plus
# one entry in `web/play.html`'s allow-list. Only these three are built for the
# web; the other examples are untouched and stay native-only.
#
# `PYTHON` exists because the renderer needs two pip packages and a PEP 668
# distribution will not install them into the system interpreter. Point it at a
# venv's python: `PYTHON=~/.venvs/isomesh/bin/python scripts/build_web.sh`.

set -euo pipefail
cd "$(dirname "$0")/.."

DEMOS=(game_mirror_dedup game_edit_tape_trim shifted_linear_root)
OUT=web/dist
PYTHON="${PYTHON:-python3}"

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

rustup target add wasm32-unknown-unknown

rm -rf "$OUT"
mkdir -p "$OUT/play/pkg"

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
echo "site total:"
du -sh "$OUT" | sed 's/^/  /'
