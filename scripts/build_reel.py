#!/usr/bin/env python3
"""Build the hosted demo reel: one self-contained HTML page, no external requests.

    scripts/build_reel.py [out.html]

`bevy_isomesh/DEMOS.md` is the in-repo demo page and renders on GitHub, crates.io
and docs.rs by pointing at the committed GIFs. This builds the *hosted* version
of the same content, which has to be one file with nothing fetched at run time --
so every clip is inlined as a `data:` URI.

# Why WebP and not the committed GIFs

The GIFs total 37 MB, and a hosted page has a 16 MB ceiling that base64 inflates
into by a third. Re-encoded to animated WebP at 620px and 14 fps the same 24
clips come to **3.6 MB**, which is not a compromise -- WebP has an alpha channel,
inter-frame prediction and no 256-colour palette, so it is both smaller and
better-looking than the GIF it came from. The GIFs stay because a GIF renders
inside a README and a WebP does not always.

Regenerate the WebP set first:

    for f in docs/gifs/*.gif; do
      ffmpeg -y -i "$f" -vf "fps=14,scale=620:-1:flags=lanczos" \
        -c:v libwebp_anim -q:v 74 -loop 0 -an "docs/gifs/web/$(basename "$f" .gif).webp"
    done
"""
import base64, pathlib, html, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
WEB = ROOT / "docs" / "gifs" / "web"
OUT = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "target" / "isomesh-reel.html"


def clip(stem):
    data = base64.b64encode((WEB / f"{stem}.webp").read_bytes()).decode()
    return f"data:image/webp;base64,{data}"

# id, tier, title, alt, prose, command
SECTIONS = [
    ("See it work", "A voxel world that has a roof, streams, and holds up under a body walking on it.", [
        ("M-69", "measured", "A world with a roof over your head", "flying-through-the-rock",
         "A camera flying through rock — under arches, into tunnels, out the far side.",
         "The field is nine lines and nothing in it is authored: <code>max(p.y − height(x,z), |gyroid(p)| − thickness)</code>. A <code>max</code> is an intersection, so rock exists only where a point is below the terrain <em>and</em> inside a thickened gyroid. A heightfield stores one number per column and cannot represent any of it.",
         "cargo run --example game_showcase --release"),
        ("M-106", "measured", "495 seam crossings, zero holes", "walking-the-seams",
         "A ball walking across streamed terrain, chunks loading continuously around it.",
         "Chunks are meshed independently, so whether two of them actually meet is decided by the overlap the layout chose. Get it wrong and a player falls through the world at a boundary — invisible in every screenshot, fatal in every playthrough. So this counts rather than asserts: <strong>495 crossings tested, 0 probes that hit nothing</strong>, worst step across a seam 0.412 cells against 0.539 within one chunk. The seams are smoother than the terrain they join.",
         "cargo run --example game_walk --release"),
        ("—", "demo", "A world larger than memory", "terrain-streaming",
         "Endless fBm terrain streaming past a flying camera.",
         "The terrain is a function, so there is no edge to reach. Fly long enough and every chunk you can see was meshed while you were watching. The number that matters is frame time <em>while chunks are landing</em>, not after.",
         "cargo run --example game_terrain_stream --release"),
        ("—", "demo", "Building the field, not meshing it", "building-a-field",
         "A capsule, a sphere, a torus — then the mushroom they add up to.",
         "The only demo with no meshing content at all: the extractor is the default, the resolution is fixed, and the one thing that changes is the expression. It shows the parts first, then <code>Union { SmoothUnion { stem, Difference { cap, flat }, k }, gills }</code>, then sweeps <code>k</code> — the stem meets the cap in a crease at zero and a fillet by 0.34. That last knob is the one a level designer reaches for, and none of this is authored geometry.",
         "cargo run --example sdf_authoring --release"),
    ]),
    ("Algorithms, side by side", "Same field, same grid, same crossings. What differs is one function.", [
        ("M-54", "measured", "The corner, and the 101×", "dual-contouring-vs-surface-nets",
         "Surface Nets rounding a box corner beside Dual Contouring holding it.",
         "The only difference between these two meshes is where a cell's vertex goes: Surface Nets takes the centroid of the crossings, Dual Contouring solves for the point that best fits the crossing planes. On a sharp field that is worth <strong>101×</strong> in Hausdorff distance — 7.217e-2 against 7.145e-4. On a smooth one it is worth 1.2×, which is why this is a choice and not an upgrade.",
         "cargo run --example dual_contouring_cube --release"),
        ("✗1", "falsified", "Fewer triangles, exactly", "surface-nets-vs-marching-cubes-box",
         "The same box under Marching Cubes and Surface Nets.",
         "Folklore says Surface Nets produces &ldquo;substantially fewer triangles&rdquo;. It produces <code>2χ</code> fewer — <code>V_sn = V_mc + χ</code> and <code>F_sn = F_mc + 2χ</code>, exactly, at every resolution. On a sphere that is four. The HUD computes it live so you can watch the identity hold.",
         "cargo run --example surface_nets_vs_marching_cubes --release"),
        ("M-67", "measured", "0 triangles against 1,340", "subgrid-letters-thinner-than-a-voxel",
         "The word ISO thinning below one voxel; Marching Cubes loses it entirely.",
         "A feature thinner than a voxel does not exist to a method that asks <em>what sign is this grid corner</em>. A sign test cannot distinguish <strong>95.6%</strong> of the configurations a tetrahedron can be in. Subgrid Marching Tetrahedra asks instead for every zero along the edge. It is also 196× the cost of Marching Cubes, which is the trade.",
         "cargo run --example subgrid_features --release"),
    ]),
    ("Game-shaped", "Chunked, edited, budgeted, collided against.", [
        ("M-33", "measured", "Digging, and what it re-meshes", "digging-a-tunnel",
         "A first-person camera carving tunnels through chunked terrain.",
         "The number this exists for is <strong>E1</strong>: a brush changes only 15–36% of the cells inside its own bounding box. Counting value changes overstates the re-mesh set by 2.8–3.7×, which is the difference between a budget that holds and one that does not.",
         "ISOMESH_AUTOCARVE=240 cargo run --example game_dig --release"),
        ("—", "demo", "The debris is the boolean", "the-debris-is-the-boolean",
         "A hollow shell being shot, cratering, debris falling away.",
         "A shot appends a subtract brush to crater the target <em>and</em> meshes the intersection of the pre-shot solid with that same sphere as the debris. Nothing is authored: the crater and the fragment are two views of one boolean.",
         "cargo run --example game_destruction --release"),
        ("M-36…38", "measured", "Undo without a snapshot", "undo-is-a-refold",
         "A solid morphing backwards and forwards as a log cursor moves.",
         "The edits are a log and the field is a fold of that log. Undo moves the cursor back and re-folds — nothing is stored, nothing copied. The order is load-bearing and the crate knows how much: same-kind hard edits are bit-identical across all 40,320 orderings, mixed add-and-subtract gives <strong>11</strong> distinct results, and smooth union gives <strong>40,317</strong>.",
         "cargo run --example game_editor --release"),
        ("—", "demo", "Paint that survives the wall", "paint-that-survives-the-wall",
         "Graffiti sprayed on a wall, then a hole blown through it.",
         "Spray, then blow a hole through the wall. The paint on what remains is exactly where you put it. The drift readout is continuously zero, because paint lives in the edit log and was never on the surface.",
         "ISOMESH_AUTOPAINT=240 cargo run --example game_paint --release"),
    ]),
    ("Where it breaks, and what that proves", "The demos that exist to fail, and the numbers that come out of them.", [
        ("✗19", "falsified", "The paper says this is a manifold", "manifold-check-resolution",
         "Non-manifold edges drawn in red on the mesh, appearing and vanishing with resolution.",
         "Defects drawn in place rather than counted in a corner. Schaefer, Ju &amp; Warren state that the uniform-grid dual <em>&ldquo;is always a manifold&rdquo;</em>. Measured over eight fields it is not — Marching Cubes gives 0 non-manifold edges and Manifold Dual Contouring gives 114. Their premise holds; their conclusion does not, and the counterexample now fits in 48 samples.",
         "cargo run --example manifold_check --release"),
        ("M-40", "measured", "The ambiguous face is rarer than you were told", "ambiguous-faces-are-rare",
         "Cells with ambiguous faces boxed in amber and magenta.",
         "Amber where the asymptotic decider agreed with plain Marching Cubes, magenta where it disagreed. Magenta is rare, and the rarity is the finding: on <strong>five of the eight</strong> reference fields the ambiguous face never occurs at all, so MC33 and plain MC are bit-identical on them at every resolution tested.",
         "cargo run --example marching_cubes_ambiguity --release"),
        ("M-222", "measured", "One cell, meshed twice", "the-tunnel-meshed-as-a-tunnel",
         "One cell as two separate discs, and as a single tunnel.",
         "Left, the face rule alone: two discs, two components, χ = 2. Right, the same cell with the interior rule: one cylinder, one component, χ = 0. A tunnel is a handle and a handle costs a closed surface exactly two — arithmetic, not opinion. The gold ring is the inner hexagon the whole construction is built from.",
         "cargo run --example marching_cubes_tunnel --release"),
        ("M-28", "measured", "What the clamp cannot reach", "qef-clamp-self-intersections",
         "Self-intersecting triangles in red, vanishing and persisting as the clamp toggles.",
         "Confining each solved vertex to its own cell drives self-intersections to <strong>exactly zero</strong> on five of the eight fields. On gyroid and fbm terrain it does not — 3.12 and 13.84 pairs per 1,000 triangles survive. What is left is a connectivity failure rather than a placement one, which is why the manifold variant makes it <em>worse</em>, not better.",
         "cargo run --example qef_clamp --release"),
        ("M-30", "measured", "The sharpness knob, at both ends", "sharp-features-lambda-sweep",
         "A box with a sphere bitten out, its edges rounding over and snapping back.",
         "λ, the Tikhonov regularizer in the vertex solve, was a compile-time constant until a demo needed to turn it. Toward zero, corners come out exactly and flat cells fling their vertices out — one landed <strong>3.18 cells</strong> outside the cell that produced it. Toward large λ, every edge rounds into Surface Nets.",
         "cargo run --example sharp_features --release"),
        ("M-32", "measured", "Where f32 blurs, and where it tears", "precision-f32-tears",
         "The same field at a large coordinate offset in f32 and f64.",
         "Two failures, two laws, and neither is where you would guess. At 10⁶ <code>f32</code> does not crack — what moves is accuracy. The crack is an order of magnitude further out, and at <strong>2²³</strong> it tears: χ drops from 2 to 1, the vertex count collapses, and real holes open.",
         "cargo run --example precision_f32_vs_f64 --release"),
        ("M-229", "measured", "The saddle, swept to infinity and back", "the-interior-decider-sweep",
         "A plane sweeping one cell while a dot traces a hyperbola.",
         "No mesh at all — one wire cell and the bilinear saddle tracked as a plane sweeps through it. The trail draws a hyperbola that runs off to infinity and returns from the other side at the pole plane. Two verdicts print side by side, and <strong>12.6%</strong> of the time the classic numerator-only test is the wrong one.",
         "cargo run --example marching_cubes_interior --release"),
    ]),
    ("On the GPU", "Faster above about 33 samples per axis, and by 37× at 129³ — provided the field is evaluated on the GPU rather than uploaded to it.", [
        ("GPU-011b", "measured", "Never touching the bus", "gpu-resident-mesh-shader",
         "A field with three brushes moving through it, extracted and drawn entirely on the GPU.",
         "Field evaluation, extraction and draw, all device-side. Per frame the CPU sends a camera matrix and three brushes, and waits. No vertex buffer, no index buffer, and zero mesh entities. <strong>0.54 ms at 129³ against a single-threaded CPU's 20.14</strong> — and that GPU figure is nearly flat across a 420× rise in cell count, because extraction was never the cost. <code>count + emit</code> is 0.045 ms and does not move. Everything the path costs is data movement, which is why where you evaluate the field decides the rest: upload CPU-sampled data instead and the upload alone is 87% of it.",
         "cargo run --example gpu_mesh_shader --release"),
        ("—", "demo", "A boolean re-meshed every frame", "a-boolean-remeshed-every-frame",
         "A CSG solid whose cutter orbits, re-meshed continuously.",
         "A concave edge is where a CAD tool lives, and it is the case where the vertex solve wants to sit outside the cell that produced it. Here it also moves. The demo reports the worst error over the whole sweep, because a single-position measurement of a sharp feature is a measurement of that position.",
         "cargo run --example game_csg_props --release"),
    ]),
]

def entry(e):
    fid, tier, title, stem, alt, prose, cmd = e
    return f"""      <article class="demo">
        <div class="rail">
          <span class="fid">{html.escape(fid)}</span>
          <span class="tier tier-{tier}">{tier}</span>
        </div>
        <div class="body">
          <h3>{title}</h3>
          <figure>
            <img src="{clip(stem)}" alt="{html.escape(alt)}" loading="lazy" decoding="async" width="620">
            <figcaption>{alt}</figcaption>
          </figure>
          <p>{prose}</p>
          <pre><code>{html.escape(cmd)}</code></pre>
        </div>
      </article>"""

sections = "\n".join(
    f"""    <section class="group">
      <header class="group-head">
        <h2>{title}</h2>
        <p>{blurb}</p>
      </header>
{chr(10).join(entry(e) for e in entries)}
    </section>"""
    for title, blurb, entries in SECTIONS)

n_clips = sum(len(entries) for _, _, entries in SECTIONS)

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(f"""<title>The isomesh Reel</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
:root {{
  --ground:#f2f3f6; --panel:#ffffff; --ink:#14161b; --muted:#5c6270;
  --rule:#dee1e8; --accent:#0e8f61; --amber:#9c6a26; --crimson:#bb3b3b;
  --shadow:0 1px 2px rgba(20,22,27,.06), 0 8px 24px rgba(20,22,27,.05);
  --mono:ui-monospace,"SF Mono","JetBrains Mono","Cascadia Code",Menlo,Consolas,monospace;
  --sans:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,"Helvetica Neue",sans-serif;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --ground:#15171c; --panel:#1c1f26; --ink:#e7e9ef; --muted:#8b91a1;
    --rule:#2a2e37; --accent:#3ddc97; --amber:#d9a45f; --crimson:#e05c5c;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 10px 30px rgba(0,0,0,.35);
  }}
}}
:root[data-theme="dark"] {{
  --ground:#15171c; --panel:#1c1f26; --ink:#e7e9ef; --muted:#8b91a1;
  --rule:#2a2e37; --accent:#3ddc97; --amber:#d9a45f; --crimson:#e05c5c;
  --shadow:0 1px 2px rgba(0,0,0,.4), 0 10px 30px rgba(0,0,0,.35);
}}

* {{ box-sizing:border-box; }}
body {{
  margin:0; background:var(--ground); color:var(--ink);
  font-family:var(--sans); font-size:16px; line-height:1.65;
  -webkit-font-smoothing:antialiased;
}}
a {{ color:var(--accent); }}
a:focus-visible, summary:focus-visible {{ outline:2px solid var(--accent); outline-offset:3px; }}

/* ---- hero: a sampling lattice, drawn rather than pictured ---------------- */
.hero {{
  position:relative; overflow:hidden;
  border-bottom:1px solid var(--rule); background:var(--panel);
}}
.hero::before {{
  content:""; position:absolute; inset:0; pointer-events:none;
  background:
    repeating-linear-gradient(0deg, var(--rule) 0 1px, transparent 1px 34px),
    repeating-linear-gradient(90deg, var(--rule) 0 1px, transparent 1px 34px);
  -webkit-mask-image:radial-gradient(120% 90% at 50% 0%, #000 10%, transparent 72%);
  mask-image:radial-gradient(120% 90% at 50% 0%, #000 10%, transparent 72%);
  opacity:.7;
}}
.hero-inner {{
  position:relative; max-width:74rem; margin:0 auto;
  padding:clamp(3rem,7vw,5.5rem) clamp(1.25rem,4vw,3rem) clamp(2.5rem,5vw,4rem);
  display:grid; gap:2.25rem;
}}
.eyebrow {{
  font-family:var(--mono); font-size:.72rem; letter-spacing:.16em;
  text-transform:uppercase; color:var(--muted); margin:0;
}}
h1 {{
  font-family:var(--mono); font-weight:700; margin:.4rem 0 0;
  font-size:clamp(2.1rem,5.5vw,3.4rem); letter-spacing:-.035em;
  line-height:1.05; text-wrap:balance;
}}
h1 .dim {{ color:var(--muted); }}
.thesis {{ max-width:34rem; margin:1rem 0 0; font-size:1.08rem; color:var(--muted); }}
.thesis strong {{ color:var(--ink); font-weight:600; }}

.figures {{
  display:grid; grid-template-columns:repeat(auto-fit,minmax(9rem,1fr));
  gap:1px; background:var(--rule); border:1px solid var(--rule); border-radius:2px;
}}
.figure {{ background:var(--panel); padding:1rem 1.1rem; }}
.figure b {{
  display:block; font-family:var(--mono); font-size:1.5rem; font-weight:700;
  letter-spacing:-.03em; font-variant-numeric:tabular-nums; color:var(--accent);
}}
.figure span {{ display:block; font-size:.82rem; color:var(--muted); margin-top:.15rem; }}

/* ---- body --------------------------------------------------------------- */
main {{ max-width:74rem; margin:0 auto; padding:0 clamp(1.25rem,4vw,3rem) 5rem; }}
.group {{ padding-top:clamp(2.75rem,6vw,4.5rem); }}
.group-head {{ border-top:2px solid var(--ink); padding-top:.9rem; margin-bottom:2rem; }}
.group-head h2 {{
  font-family:var(--mono); font-size:1.22rem; font-weight:700;
  letter-spacing:-.02em; margin:0;
}}
.group-head p {{ margin:.3rem 0 0; color:var(--muted); font-size:.95rem; }}

.demo {{
  display:grid; grid-template-columns:7.5rem minmax(0,1fr);
  gap:0 2rem; padding:1.9rem 0; border-top:1px solid var(--rule);
}}
.demo:first-of-type {{ border-top:none; padding-top:.25rem; }}
.rail {{ display:flex; flex-direction:column; gap:.4rem; align-items:flex-start; }}
.fid {{
  font-family:var(--mono); font-size:.85rem; font-weight:700;
  letter-spacing:-.01em; font-variant-numeric:tabular-nums;
}}
.tier {{
  font-family:var(--mono); font-size:.62rem; letter-spacing:.13em;
  text-transform:uppercase; padding:.16rem .45rem; border:1px solid currentColor;
  border-radius:2px; line-height:1.5;
}}
.tier-measured {{ color:var(--accent); }}
.tier-falsified {{ color:var(--crimson); }}
.tier-demo {{ color:var(--muted); }}

.body h3 {{
  font-family:var(--mono); font-size:1.02rem; font-weight:700; margin:0 0 .85rem;
  letter-spacing:-.015em; text-wrap:balance;
}}
.body p {{ margin:1rem 0 0; max-width:44rem; }}
figure {{ margin:0; }}
figure img {{
  display:block; width:100%; max-width:620px; height:auto;
  border:1px solid var(--rule); border-radius:2px; box-shadow:var(--shadow);
  background:#1e2128;
}}
figcaption {{
  font-size:.82rem; color:var(--muted); margin-top:.5rem;
  max-width:44rem; font-style:italic;
}}
pre {{
  margin:1rem 0 0; padding:.7rem .9rem; overflow-x:auto;
  background:var(--ground); border:1px solid var(--rule); border-radius:2px;
}}
pre code {{ font-family:var(--mono); font-size:.82rem; color:var(--ink); }}
code {{ font-family:var(--mono); font-size:.9em; }}
p code {{ color:var(--amber); }}

footer {{
  border-top:1px solid var(--rule); background:var(--panel);
  padding:2.5rem clamp(1.25rem,4vw,3rem);
}}
.foot-inner {{
  max-width:74rem; margin:0 auto; display:flex; flex-wrap:wrap;
  gap:1.5rem 2.5rem; align-items:baseline;
}}
.foot-inner p {{ margin:0; color:var(--muted); font-size:.9rem; max-width:32rem; }}
.foot-links {{ display:flex; flex-wrap:wrap; gap:1.25rem; font-family:var(--mono); font-size:.85rem; }}

@media (max-width:820px) {{
  .demo {{ grid-template-columns:1fr; gap:.75rem; }}
  .rail {{ flex-direction:row; align-items:center; gap:.6rem; }}
}}
@media (prefers-reduced-motion: reduce) {{
  * {{ animation:none !important; transition:none !important; }}
}}
</style>

<div class="hero">
  <div class="hero-inner">
    <div>
      <p class="eyebrow">isomesh · engine-agnostic isosurface extraction in Rust</p>
      <h1>Twenty-four demos,<br><span class="dim">each with its receipt.</span></h1>
      <p class="thesis">Every clip below is <code>cargo run --example</code> on a machine you can
      reproduce. Every claim beside one carries the finding that produced it — <strong>measured
      here</strong>, or <strong>falsified</strong> and kept anyway, because which sources to distrust
      is worth more than the individual fact.</p>
    </div>
    <figure>
      <img src="{clip('kitchen-sink')}" alt="Eight isomesh examples running at once in a four-by-two grid" width="620">
      <figcaption>Eight of them at once. Top: a field assembled from primitives; caves and arches from
      nine lines; letters thinner than one voxel; one cell meshed as two discs and as a tunnel.
      Bottom: Surface Nets against Marching Cubes; a tunnel being carved; 495 seam crossings walked;
      a slab blown apart, where the debris is the boolean.</figcaption>
    </figure>
    <div class="figures">
      <div class="figure"><b>4.26&times;</b><span>faster dual meshing, byte-identical output</span></div>
      <div class="figure"><b>363</b><span>entries in the findings ledger</span></div>
      <div class="figure"><b>19</b><span>beliefs falsified, none deleted</span></div>
      <div class="figure"><b>17</b><span>predictions registered before measuring</span></div>
      <div class="figure"><b>1</b><span>dependency in the core crate</span></div>
    </div>
  </div>
</div>

<main>
{sections}
</main>

<footer>
  <div class="foot-inner">
    <p><strong>Vibe coded.</strong> Written by an AI agent working from a research corpus and a ticket
    queue. Every number here is produced by a test in the repository and every algorithm cites its
    source, but none of it has been through human code review. Read the tests first.</p>
    <nav class="foot-links">
      <a href="https://github.com/ladvien/isomesh">github</a>
      <a href="https://crates.io/crates/isomesh">crates.io</a>
      <a href="https://docs.rs/isomesh">docs.rs</a>
      <a href="https://github.com/ladvien/isomesh/blob/main/docs/experiments.md">the experiments</a>
      <a href="https://github.com/ladvien/isomesh/blob/main/FINDINGS.md">the ledger</a>
    </nav>
  </div>
</footer>
""")
size = OUT.stat().st_size
print(f"wrote {OUT} — {size/1048576:.2f} MB, {n_clips} demo clips + hero")
