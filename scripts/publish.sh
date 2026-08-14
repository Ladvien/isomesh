#!/usr/bin/env bash
# Publish every workspace crate whose manifest version is not yet on crates.io.
#
#   ./scripts/publish.sh --dry-run   # package and verify, upload nothing
#   ./scripts/publish.sh             # upload what is missing
#
# The logic lives here rather than in the workflow for the same reason
# `backlog_gate.sh` does: a developer and CI run the same thing, so "it passed
# locally" and "it passed in CI" mean the same sentence.
#
# # Why this is version-driven rather than push-driven
#
# Publishing cannot be undone. A version can be yanked, never deleted, and the
# name is claimed forever. So this does not publish because a push happened; it
# publishes because a version exists in the manifest that does not exist on the
# registry. Pushing a version bump releases. Pushing anything else is a no-op
# that exits 0 — the alternative, failing on "crate version already exists",
# would leave main permanently red and train everyone to ignore it.
#
# # Why the order is written down rather than derived
#
# `isomesh-gpu` depends on `isomesh` by path *and version*, and that version is
# resolved against crates.io during its own verification build. Publish it first
# and it fails with "failed to select a version for the requirement
# isomesh = ^x.y.z" — measured, not predicted. Two crates do not justify a
# topological sort, but they do justify the guard below: a new workspace member
# that nobody placed in this list is an error, not a silent omission.
set -euo pipefail

# Dependency order. Every workspace member must appear.
ORDER=(isomesh isomesh-gpu)

DRY_RUN=0
if [ "${1:-}" = "--dry-run" ]; then
  DRY_RUN=1
elif [ -n "${1:-}" ]; then
  echo "usage: $0 [--dry-run]" >&2
  exit 2
fi

cd "$(dirname "$0")/.."

# Every crate cargo believes is in this workspace, so the list above cannot
# quietly fall behind the manifest.
#
# `--no-deps` is what makes `.packages` mean "workspace members" rather than the
# whole resolved graph, and jq is what makes it correct: the first version of
# this grepped `"name":"..."` out of the raw JSON and reported `criterion` as a
# workspace member, because a dev-dependency's name is also a `"name"` field.
mapfile -t MEMBERS < <(
  cargo metadata --format-version 1 --no-deps | jq -r '.packages[].name' | sort -u
)
for member in "${MEMBERS[@]}"; do
  found=0
  for named in "${ORDER[@]}"; do
    [ "$member" = "$named" ] && found=1
  done
  if [ "$found" -eq 0 ]; then
    echo "::error::workspace member '$member' is not in publish.sh's ORDER list."
    echo "Add it in dependency order -- publishing a crate before what it depends on fails."
    exit 1
  fi
done

published=0
skipped=0

for crate in "${ORDER[@]}"; do
  version=$(
    cargo metadata --format-version 1 --no-deps \
      | jq -r --arg c "$crate" '.packages[] | select(.name == $c) | .version'
  )
  if [ -z "$version" ]; then
    echo "::error::could not read a version for $crate"
    exit 1
  fi

  # 404 means "not published". Anything else -- including a network failure --
  # must not be read as "not published", because that would republish under a
  # version that already exists and fail confusingly, or worse, succeed against
  # the wrong registry.
  code=$(curl -sS -o /dev/null -w '%{http_code}' \
    -H 'User-Agent: isomesh-publish (github actions)' \
    "https://crates.io/api/v1/crates/$crate/$version")

  case "$code" in
    200)
      echo "== $crate $version is already on crates.io -- nothing to do"
      skipped=$((skipped + 1))
      continue
      ;;
    404)
      echo "== $crate $version is not on crates.io"
      ;;
    *)
      echo "::error::crates.io returned HTTP $code for $crate $version; refusing to guess"
      exit 1
      ;;
  esac

  if [ "$DRY_RUN" -eq 1 ]; then
    # `--dry-run` packages and verification-builds without uploading. It is the
    # step that catches the errors which are only visible at package time --
    # a `path` dependency with no `version`, a missing license file, a `README`
    # outside the package directory -- none of which any test can see.
    #
    # It cannot verify `isomesh-gpu` before `isomesh` is on the registry, so on
    # a first release of a dependent crate this reports that and moves on
    # rather than failing the run.
    echo "-- dry run: cargo publish -p $crate --dry-run"
    if ! cargo publish -p "$crate" --dry-run --locked; then
      echo "::warning::$crate could not be verified. If it depends on another"
      echo "crate in this workspace whose version is not yet on crates.io, that"
      echo "is expected and resolves once the dependency is published."
    fi
    continue
  fi

  # **The token is required here, and deliberately not a line earlier.** An
  # upload is now imminent -- this crate's version is absent from the registry
  # and this is not a dry run -- so a missing token means a release that will
  # silently not happen, and it must fail loudly.
  #
  # It used to be checked at the top of the workflow job instead, which made
  # *every* push to a green main red, including the ones that upload nothing.
  # That is precisely the outcome this script's own header exists to prevent:
  # "would leave main permanently red and train everyone to ignore it". Nobody
  # had seen it because the job had been `skipped` on every run since it was
  # written -- the suite was red -- so GPU-013's push was the first time it ran
  # at all, and the first time the missing secret was visible (M-198).
  if [ -z "${CARGO_REGISTRY_TOKEN:-}" ]; then
    echo "::error::CARGO_REGISTRY_TOKEN is not set, and $crate $version needs uploading."
    echo "Add it under Settings -> Environments -> crates-io -> Environment secrets."
    exit 1
  fi

  echo "-- cargo publish -p $crate"
  # `--locked` so CI publishes exactly the dependency graph the committed
  # lockfile describes, rather than whatever resolved that morning.
  cargo publish -p "$crate" --locked
  published=$((published + 1))
done

echo
if [ "$DRY_RUN" -eq 1 ]; then
  echo "dry run complete: $skipped already published, $((${#ORDER[@]} - skipped)) would be uploaded"
else
  echo "published $published, skipped $skipped already on crates.io"
fi
