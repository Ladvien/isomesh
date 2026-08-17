#!/usr/bin/env bash
#
# Fetch the Open SciVis volumes the comparison benchmarks read.
#
# Every reference field in this crate is **analytic**, which is what makes the
# accuracy harness exact and is exactly why its timings cannot be set beside a
# published isosurfacing benchmark -- those are run on real CT and simulation
# volumes. M-006. These are those volumes.
#
# The data is NOT committed. It is downloaded here, into a gitignored directory,
# and verified against **the publisher's own SHA-512**, which the dataset page
# lists per file. That distinction matters: a hash this script computed from the
# file it had just downloaded would verify nothing at all.
# `docs/measurements/volumes/PROVENANCE.md` is committed
# and carries the URL, hash and licence of each file, so the provenance survives
# in git while the bytes do not.
#
# ## The transport is HTTP, and that is the server's choice rather than ours
#
# klacansky.com serves Open SciVis Datasets over **HTTP only** -- port 443
# refuses the connection outright (measured, V-40), so there is no HTTPS URL to
# prefer and nothing to fall back from. Integrity comes from the published
# SHA-512 below, not from the transport, which is the right guarantee for
# content-addressed data anyway: a wrong or tampered file fails the hash, is
# deleted, and this script exits non-zero.
#
# ## Usage
#
#   ./scripts/fetch_volumes.sh          # fetch anything missing, verify all
#   ./scripts/fetch_volumes.sh --check  # verify what is present, fetch nothing
#
# Benchmarks that read these skip cleanly when they are absent, so a clean clone
# with no network still builds and tests.

set -euo pipefail

cd "$(dirname "$0")/.."

DEST=docs/measurements/volumes
BASE=http://klacansky.com/open-scivis-datasets

# name | subdirectory | sha512, as published on the dataset page
#
# Two volumes, both `uint8`, chosen for what they buy rather than for fame.
#
#   fuel   -- 64^3, 256 KB. Small enough to fetch in a second and to run every
#             extractor over at full resolution, including subgrid Marching
#             Tetrahedra at ~200x Marching Cubes (M-308).
#   bonsai -- 256^3, 16 MB. The volume published comparisons actually use, which
#             is the entire point of M-006.
#
# `uint8` is deliberate on both. Quantised data is what makes Grosso 2017's
# singular faces reachable -- 8, 58 and 20 per CT volume, where a continuous f64
# field produces zero (M-220, M-232) -- so these are also the fixture A-002i and
# A-020b have been waiting for.
VOLUMES=(
    "fuel_64x64x64_uint8.raw|fuel|77fdd7c657da1946bafc84e88c6b8a03ae104a79a5bdec3c7db9257480ef4bf72551a08d22fd237c8e387dd2571b575f1a1a11f5f32b1fa4d4ef385d9fe1d613"
    "bonsai_256x256x256_uint8.raw|bonsai|b34156a0ffc80ffaf84d069f3d05a40fdd999a35f05492829a2b0c13403a3147e73712b1d10c2cc34da66a59540a1632dae6adc96f3ebf3efa5d4d6c10598997"
)

check_only=0
if [ "${1:-}" = "--check" ]; then
    check_only=1
elif [ "$#" -gt 0 ]; then
    printf 'usage: %s [--check]\n' "$0" >&2
    exit 2
fi

mkdir -p "$DEST"

fail=0
for entry in "${VOLUMES[@]}"; do
    IFS='|' read -r name subdir want <<<"$entry"
    path="$DEST/$name"

    if [ ! -f "$path" ]; then
        if [ "$check_only" -eq 1 ]; then
            printf '  missing  %s\n' "$name"
            continue
        fi
        printf '  fetching %s ... ' "$name"
        # To a temporary name, so an interrupted download never looks complete.
        if ! curl -fsS --max-time 900 -o "$path.part" "$BASE/$subdir/$name"; then
            printf 'FAILED\n'
            rm -f "$path.part"
            fail=1
            continue
        fi
        mv "$path.part" "$path"
        printf 'done\n'
    fi

    got=$(sha512sum "$path" | cut -d' ' -f1)
    if [ "$got" != "$want" ]; then
        printf '  %s\n    ::error:: sha512 mismatch against the published hash\n      expected %s\n      got      %s\n' \
            "$name" "$want" "$got"
        rm -f "$path"
        fail=1
    else
        printf '  ok       %s\n' "$name"
    fi
done

if [ "$fail" -ne 0 ]; then
    printf '\nfetch_volumes FAILED\n' >&2
    exit 1
fi

printf '\nvolumes in %s\n' "$DEST"
