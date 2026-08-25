---
name: isomesh-site-deploy
description: Ship an isomesh change to isomesh.ladvien.com and prove it landed — the pre-push gates, the push-to-main deploy, and the artefact-level verification that catches a stale wasm module behind fresh page prose. Use whenever a commit has to become live, or when asked "is it live yet".
---

# Getting an isomesh change live on isomesh.ladvien.com

There is **one** deploy path and it is a push to `main`. No manual upload, no
`gh workflow run`, no branch deploy. `.github/workflows/ci.yml`'s `pages` job is
gated `if: github.event_name == 'push' && github.ref == 'refs/heads/main'` and
`needs: site`, so nothing reaches the CDN except a green `site` build of a commit
on `main`.

Budget **~30 minutes** of CI. The `site` job installs `wasm-bindgen-cli` from the
lockfile version and builds nine Bevy wasm demos; the last measured run was
**28m32s** wall from push to `deploy to github pages` finishing.

## Before you push

Run these from the repo root. All three are the same commands CI runs, which is
the point — a local green is evidence, not a rehearsal.

```bash
./scripts/backlog_gate.sh                                    # counts vs rows, both files
PYTHON=~/.venvs/isomesh/bin/python ./scripts/preflight.sh --full
```

`preflight --full` is **4m37s** on the 5900X. Run it *supervised*, not as a
foreground `bash` call piped into `grep`: the pipe buffers the step headers, the
call blocks past the harness's foreground window, and you cannot see which step
it is on. Either let it write a log and read that, or start it as a process:

```bash
PYTHON=~/.venvs/isomesh/bin/python ./scripts/preflight.sh --full 2>&1 | tee /tmp/preflight.log
grep -aE "^── |FAILED|all green" /tmp/preflight.log
```

Expect **23 `ok` steps** and `preflight --full: all green`. Note what it does
*not* cover: `bevy: test` is `cargo test --lib`, so the **example** tests — every
`game_dig` GPU gate among them — are yours to run:

```bash
cd bevy_isomesh && cargo test --example game_dig
```

Toolchain: stable **1.98 or newer**. Clippy 0.1.98 enforces lints that do not
exist in 1.97, so a green 1.97 run is not evidence (M-304).

## Committing

`CLAUDE.md` rule 2: **one ticket, one commit, message starts with the ticket ID.**
Move the row from `BACKLOG.md` to `BACKLOG_ARCHIVE.md` in that same commit and fix
all three counts — the archive's `N tickets.`, and `BACKLOG.md`'s
`**N tickets archived, M open.**`. `backlog_gate.sh` checks them against the rows
and fails loudly on drift.

If the commit or the archive row quotes a measured number, **re-read it from the
committed tree**. A figure carried over from a build you have since changed
describes a run that no longer exists — this repository has published two tables
that cited runs absent from the files they referenced.

The same push also runs `publish to crates.io`. That is not a surprise and not a
release: `publish.sh` uploads a crate only when its manifest version is absent
from crates.io, so an ordinary push exits 0 having uploaded nothing. Only a
version bump releases.

## Push and watch

```bash
git push origin main
gh run list --branch main --limit 3 \
  --json databaseId,status,conclusion,headSha,displayTitle \
  --jq '.[] | "\(.databaseId) \(.status) \(.conclusion // "-") \(.headSha[0:7]) \(.displayTitle[0:50])"'
```

`gh run watch <id> --interval 30 --exit-status` blocks for the whole build, so run
it as a supervised process rather than a foreground command. Then confirm every
job, and specifically that the deploy ran:

```bash
gh run view <id> --json jobs --jq '.jobs[] | "\(.conclusion // .status)\t\(.name)"'
```

`deploy to github pages` must read `success`. If it says `Not Found`, the Pages
source is not set to "GitHub Actions" under Settings → Pages, and no workflow can
set it for you.

## Verifying it is live — the artefact, not the page

**Use `http://`, never `https://`.** The custom domain has no certificate:
`gh api repos/Ladvien/isomesh/pages` reports `https_certificate` **absent** and
`status: null`, so `https://isomesh.ladvien.com` fails the TLS handshake outright
while `http://` answers 200. Check before assuming it is fixed:

```bash
gh api repos/Ladvien/isomesh/pages \
  --jq '{status, https_enforced, cert: (.https_certificate.state // "absent")}'
```

Once that cert appears, one command finishes the job —
`gh api -X PUT repos/Ladvien/isomesh/pages -F https_enforced=true` — and it must
not be run before, because enforcing https without a certificate takes the site
down rather than securing it.

### 1. The rendered page

```bash
curl -sS -H 'Cache-Control: no-cache' -o /tmp/live.html \
  "http://isomesh.ladvien.com/play.html?cb=$(date +%s)"
grep -c "<the new prose you added>" /tmp/live.html
```

### 2. The compiled module — this is the check that matters

A page-prose check passes happily while the wasm behind it is stale, because the
markdown and the module are built by the same job but from different inputs. Grep
the deployed binary for a string only the new code contains **and** for the
absence of the string it replaced. Both directions, or it proves nothing:

```bash
curl -sS -o /tmp/m.wasm http://isomesh.ladvien.com/play/pkg/game_dig/game_dig_bg.wasm
strings -a /tmp/m.wasm | grep -c '<new format-string fragment>'   # expect >= 1
strings -a /tmp/m.wasm | grep -c '<old format-string fragment>'   # expect 0
```

Rust splits a format literal around its `{}` holes, so grep the **fragments**
(`" (cached)"`, `"KB of samples uploaded since start"`), never the whole
`"{:.2} ms"` template — that never appears contiguously and its absence means
nothing. Demo modules live at `web/dist/play/pkg/<demo>/<demo>_bg.wasm`; the
deployed URL mirrors that path.

## What cannot be verified on this host, and what to do instead

`big` has **no reachable display** (`XDG_SESSION_TYPE=tty`, no wayland socket),
and headless Chrome here exposes **no `navigator.gpu`** — `--enable-unsafe-webgpu`
with SwiftShader through ANGLE does not change that. Every playable demo needs
WebGPU, so the live canvas cannot be driven from this machine. What you get is the
page's own honest gate: *"This demo needs WebGPU, which this browser has not
enabled."* Seeing that text is a real result — it proves the guard ships — but it
is not the demo running.

So verify demo *behaviour* natively, against real hardware, and say plainly that
the live canvas was not driven. A HUD string is testable without a window: run
`report` as a one-shot system over the demo's own harness and read `DemoStats`
(`game_dig`'s `the_hud_reports_the_cached_base_and_the_measured_verdict` is the
pattern — it prints the whole panel and asserts the lines a reader is meant to
read).

If you do drive Chrome here, two flags are mandatory or the navigation looks
blocked when it is really the certificate:

```bash
chrome --headless=new --no-sandbox \
  --ignore-certificate-errors \
  --disable-features=HttpsUpgrades,HttpsFirstBalancedModeAutoEnable,HttpsFirstModeV2 \
  --remote-debugging-port=45123 --user-data-dir=/tmp/cdp
```

Without them Chrome silently upgrades `http://` to `https://`, the handshake fails
(`net_error -200`), and the tab reports **`ERR_BLOCKED_BY_CLIENT`** — which reads
like an extension or a policy block and is neither. Spawn Chrome yourself with an
explicit `--headless=new` and `--remote-debugging-port`, then attach the browser
tool by `cdp_url`: letting the tool spawn its own browser times out waiting for a
CDP endpoint on a display-less host.

## Known-benign findings

- `GET /favicon.ico` → **404** on every page. Cosmetic, pre-existing, not yours.
- Bevy 0.19's slab allocator logs *"Use-after-free: attempted to copy element data
  for an unallocated key"* thousands of times in `game_showcase`. Upstream,
  documented on `play.html` rather than silenced.
