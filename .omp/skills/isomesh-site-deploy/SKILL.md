---
name: isomesh-site-deploy
description: Ship an isomesh change to ladvien.github.io/isomesh and prove it landed — the pre-push gates, the push-to-main deploy, and the artefact-level verification that catches a stale wasm module behind fresh page prose. Use whenever a commit has to become live, or when asked "is it live yet".
---

# Getting an isomesh change live on ladvien.github.io/isomesh

There is **one** deploy path and it is a push to `main`. No manual upload, no
`gh workflow run`, no branch deploy. `.github/workflows/ci.yml`'s `pages` job is
gated `if: github.event_name == 'push' && github.ref == 'refs/heads/main'` and
`needs: site`, so nothing reaches the CDN except a green `site` build of a commit
on `main`.

Budget **~30 minutes** of CI. The `site` job installs `wasm-bindgen-cli` from the
lockfile version and builds nine Bevy wasm demos; the last measured run was
**28m32s** wall from push to `deploy to github pages` finishing.

## Before you push

Run these from the repo root. Both are commands CI runs too, which is the point —
a local green is evidence, not a rehearsal.

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
the counts with it: **three numbers across two files** — the archive's
`N tickets.` above its index, and both halves of `BACKLOG.md`'s
`**N tickets archived, M open.**`. `backlog_gate.sh` checks all three against the
rows and fails loudly on drift.

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

**Use `https://ladvien.github.io/isomesh/`. There is no custom domain, and that
is deliberate.** `isomesh.ladvien.com` was the custom domain and GitHub never
issued its certificate: `pages/health` reported `is_valid: true`,
`is_https_eligible: true`, `caa_error: null`, `is_proxied: false` and
`https_error: "peer_failed_verification"` while `https_certificate` stayed
**absent**, for over a day. DNS was verified independently and was clean — CNAME
to `ladvien.github.io`, GitHub's four Pages IPs, Route 53, no proxy, and
`github.io`'s CAA permits `letsencrypt.org` — so provisioning was stalled on
GitHub's side with nothing to fix locally.

That was not cosmetic. **WebGPU is secure-context-only**, every demo on the site
needs it, and `https://ladvien.github.io/isomesh/` 301-redirected to the
`http://` custom domain — so `navigator.gpu` was `undefined` and `play.html`'s
gate fired in *every* browser, desktop Chrome included. The pretty URL cost the
whole site.

**A second, non-obvious blocker, and the one that actually took the site down.**
Dropping isomesh's own custom domain was not enough: `Ladvien/ladvien.github.io`
— an archived 2022 fork — still carried Pages custom domain `ladvien.com`, and an
**account-level** custom domain makes GitHub 301 every *project* page
(`ladvien.github.io/<repo>/` → `ladvien.com/<repo>/`). `ladvien.com` is served by
AWS S3 + CloudFront (`server: AmazonS3`), never by GitHub Pages, so those
redirects landed on a 404. `gh api -X PUT` on an archived repo returns **409
"Repository is archived."**, so clearing it took three calls:

```bash
gh api -X PATCH repos/Ladvien/ladvien.github.io -F archived=false
printf '%s' '{"cname": null}' | gh api -X PUT repos/Ladvien/ladvien.github.io/pages --input -
gh api -X PATCH repos/Ladvien/ladvien.github.io -F archived=true
```

`ladvien.com` was unaffected, because CloudFront serves it. If a project page
ever 301s somewhere unexpected again, check the *user* site's `cname` first — the
repo you are deploying is not where the redirect comes from.

**`isomesh.ladvien.com` is gone from DNS too, and that was a second step.**
Clearing the Pages setting left the Route 53 `CNAME` pointing at
`ladvien.github.io`, which is the worst of the three possible states: the name
still resolved, so `https://` failed the TLS handshake (`curl` exit 60, no
certificate for that hostname) and `http://` returned GitHub's own 404 for a
`Host` it no longer maps. **That reads exactly like an outage**, and it cost a
reader a puzzled minute — the fix for which is not documentation but deleting the
record. Done in the Route 53 console, because it is a one-off deletion of one
record and the console shows you the row before you commit it; the CLI route
needs a valid key, and the one on `big` has been invalid since April.

So the expected state of the old hostname is **`NXDOMAIN`**, and that is the
check:

```bash
dig +short isomesh.ladvien.com @1.1.1.1     # expect no output
dig +short ladvien.com          @1.1.1.1     # expect CloudFront, untouched
```

An answer from the first means the record is back and the *next* deploy could
re-provision the custom domain. Note the second line: the record lived in the
`ladvien.com` zone beside the apex that serves the personal site, so "did I
delete only the one row" is worth one command.

The site now serves under GitHub's own `*.github.io` certificate, which is valid,
and `https_enforced` is `true`. Confirm both, and that no domain has crept back:

```bash
gh api repos/Ladvien/isomesh/pages --jq '{cname, html_url, https_enforced, status}'
curl -sS -o /dev/null -w '%{http_code} %{url_effective}\n' \
  https://ladvien.github.io/isomesh/play.html
```

`cname` must be `null` and the curl must print `200` and that same URL. **A `301`
means a custom domain is back** — most likely because someone re-added the
`CNAME` write to `scripts/build_web.sh`. A `CNAME` file in the artifact *re-sets*
the custom domain on every deploy, which is why that line is gone and must stay
gone.

### 1. The rendered page

```bash
curl -sS -H 'Cache-Control: no-cache' -o /tmp/live.html \
  "https://ladvien.github.io/isomesh/play.html?cb=$(date +%s)"
grep -c "<the new prose you added>" /tmp/live.html
```

### 2. The compiled module — this is the check that matters

The page prose and the demo behaviour have **different authors**: a human writes
`web/play.html`, the compiler writes the module. So a prose grep confirms only
that someone described the change — it passes just as happily over a module that
does not do it, and over a partial or cached deploy. Grep the deployed binary for
a string only the new code contains **and** for the absence of the string it
replaced. Both directions, or it proves nothing:

```bash
curl -sS -o /tmp/m.wasm https://ladvien.github.io/isomesh/play/pkg/game_dig/game_dig_bg.wasm
strings -a /tmp/m.wasm | grep -c '<new format-string fragment>'   # expect >= 1
strings -a /tmp/m.wasm | grep -c '<old format-string fragment>'   # expect 0
```

`grep -c` **exits 1 when the count is zero**, so the second command "fails" on
success. Under `set -e` that aborts the script; read the number, not the exit
status.

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
page's own honest gate — *"This demo needs WebGPU. Safari needs 26 or newer…"* if
`navigator.gpu` is missing, or *"This browser has WebGPU but no usable GPU
adapter"* if `requestAdapter()` returns null. Seeing either text is a real result
— it proves the guard ships — but it is not the demo running.

So verify demo *behaviour* natively, against real hardware, and say plainly that
the live canvas was not driven. A HUD string is testable without a window: run
`report` as a one-shot system over the demo's own harness and read `DemoStats`
(`game_dig`'s `the_hud_reports_the_cached_base_and_the_measured_verdict` is the
pattern — it prints the whole panel and asserts the lines a reader is meant to
read).

If you do drive Chrome here, disable the cache before reloading or a stylesheet
or page change silently measures the old bytes: send `Network.setCacheDisabled`
with `true` on the CDP session first.

```bash
chrome --headless=new --no-sandbox \
  --remote-debugging-port=45123 --user-data-dir=/tmp/cdp
```

Spawn Chrome yourself with an explicit `--headless=new` and
`--remote-debugging-port`, then attach the browser tool by `cdp_url`: letting the
tool spawn its own browser times out waiting for a CDP endpoint on a display-less
host. The `--ignore-certificate-errors` and `--disable-features=HttpsUpgrades,…`
flags this section used to demand are **no longer needed** — they existed only
because the site was `http://` and Chrome kept upgrading it to a hostname with no
certificate, which surfaced as a misleading `ERR_BLOCKED_BY_CLIENT`.

## Known-benign findings

- `GET /favicon.ico` → **404** on every page. Cosmetic, pre-existing, not yours.
- Bevy 0.19's slab allocator logs *"Use-after-free: attempted to copy element data
  for an unallocated key"* thousands of times in `game_showcase`. Upstream,
  documented on `play.html` rather than silenced.
