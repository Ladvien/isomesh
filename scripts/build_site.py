"""Render isomesh's public prose into `web/dist`, and copy its assets beside it.

Run as `python3 scripts/build_site.py`, from anywhere: paths resolve against this
file rather than the working directory. It creates and overwrites but never
deletes, so it can be re-run on its own while iterating on prose without
rebuilding any of the three ~30 MB wasm demos. `scripts/build_web.sh` owns the
clean, and is the entry point that builds both halves.

# It is also the link checker

Every relative link and image target in the seven rendered sources is resolved
against the repository, and a target that resolves to nothing fails the build with
the source that carries it. That is the whole reason rendering only part of the
repository is safe: `BACKLOG.md`, `docs/research/*` and `crates/isomesh/src/*`
are not rendered, so links into them are rewritten to github.com blob URLs -- and
if one of those paths is ever deleted or moved, this build says so rather than
publishing a dead link.
"""

import re
import shutil
import sys
from pathlib import Path
from typing import NamedTuple

try:
    import markdown
    import pymdownx  # noqa: F401  -- presence check for the tilde extension
except ImportError:
    sys.exit(
        "build_site.py needs two packages:\n"
        "  python3 -m pip install --user 'markdown>=3.7,<4' "
        "'pymdown-extensions>=10,<11'\n"
        "On a PEP 668 distribution that pip refuses; make a venv and use its\n"
        "python, or set PYTHON=/path/to/venv/bin/python for build_web.sh."
    )

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "web" / "dist"

# What the site contains, in nav order. `nav_label` of None keeps a page out of
# the nav bar but still renders it and still lets other pages link to it.
#
# `link_base` is the directory a page's relative links resolve against, and it is
# the column an implementer is most likely to assume away: `docs/demos/*.md`
# spell their screenshots `../screenshots/…`, 26 distinct files of them, and
# treating those as repo-root-relative breaks every one silently.
PAGES = [
    ("web/index.md", "index.html", "Home", ""),
    ("FINDINGS.md", "findings.html", "Findings", ""),
    ("bevy_isomesh/DEMOS.md", "demos.html", "Demos", "bevy_isomesh"),
    ("docs/demos/algorithms.md", "demos-algorithms.html", None, "docs/demos"),
    ("docs/demos/correctness.md", "demos-correctness.html", None, "docs/demos"),
    ("docs/demos/gameplay.md", "demos-gameplay.html", None, "docs/demos"),
    ("docs/experiments.md", "experiments.html", "Experiments", "docs"),
    ("README.md", "readme.html", "Readme", ""),
]

# Copied wholesale rather than filtered to the referenced files. Every GIF is
# referenced from `DEMOS.md` anyway, and for screenshots the gap is 26 referenced
# of 51 -- 18 MB of a 200 MB site, which is not worth a filter that can be wrong.
ASSET_DIRS = [
    ("docs/gifs", "gifs"),
    ("docs/screenshots", "screenshots"),
    ("docs/experiments", "experiments"),
]

REPO = "https://github.com/ladvien/isomesh"
RAW = "https://raw.githubusercontent.com/ladvien/isomesh/main/"

# `extra` supplies tables, fenced code, attribute lists, footnotes and definition
# lists; `toc` supplies heading ids, which one link in `docs/experiments.md`
# depends on. `pymdownx.tilde` is why the second package is a dependency:
# Python-Markdown has no strikethrough and `FINDINGS.md` uses `~~…~~` fifteen
# times for genuine prose strikethrough -- retracted claims, falsified bounds.
EXTENSIONS = ["extra", "sane_lists", "toc", "pymdownx.tilde"]

BY_SOURCE = {source: output for source, output, _, _ in PAGES}

# Collected rather than raised on, so one run names every broken link instead of
# the first one.
FAILURES: list[tuple[str, str]] = []


def nav() -> str:
    """The nav bar, from the `nav_label` column plus a link to the repository."""
    links = [
        f'<a href="{output}">{label}</a>'
        for _, output, label, _ in PAGES
        if label is not None
    ]
    links.append(f'<a href="{REPO}">GitHub</a>')
    return "\n".join(links)


def rewrite(target: str, source: str, link_base: str) -> str:
    """Point one markdown link or image target at its place in the site.

    The rules, in this order and no other:

    1. A `#fragment` is split off before anything else and re-attached at the
       end. Two targets in the corpus carry one and would otherwise be resolved
       as a path that does not exist.
    2. `site:` means a site path with no repository counterpart -- how
       `web/index.md` links `play.html?demo=…`. The remainder is emitted
       verbatim.
    3. A bare fragment is left alone; `toc` has produced the ids.
    4. A `raw.githubusercontent.com` URL is stripped back to a repository path
       and falls through to rule 6. This is what localises the 55 hotlinked GIFs.
    5. Any other absolute URL is left alone.
    6. A repository path that is a rendered source becomes that page.
    7. A repository path under an asset directory becomes the copied asset.
    8. Any other path that exists becomes a github.com blob or tree URL.
    9. Anything left is a broken link, and is collected.
    """
    path, _, fragment = target.partition("#")
    fragment = f"#{fragment}" if fragment else ""

    if target.startswith("site:"):
        return target.removeprefix("site:")
    if not path:
        return target
    if path.startswith(RAW):
        # Repo-root-relative, not `link_base`-relative: the prefix already spelled
        # the whole path from the root. `DEMOS.md` hotlinks 34 GIFs and resolves
        # its own relative links against `bevy_isomesh/`, so conflating the two
        # turns every hotlink into `bevy_isomesh/docs/gifs/…` and loses all 55.
        path = path.removeprefix(RAW)
        link_base = ""
    elif path.startswith(("http://", "https://", "mailto:")):
        return target

    resolved = str((Path(link_base) / path).as_posix())
    # `Path` cannot normalise `..` on its own without touching the filesystem.
    parts: list[str] = []
    for part in resolved.split("/"):
        if part in ("", "."):
            continue
        if part == ".." and parts and parts[-1] != "..":
            parts.pop()
        else:
            parts.append(part)
    resolved = "/".join(parts)

    if resolved in BY_SOURCE:
        return BY_SOURCE[resolved] + fragment
    for source_dir, dest in ASSET_DIRS:
        if resolved.startswith(f"{source_dir}/"):
            return f"{dest}/{resolved.removeprefix(source_dir + '/')}{fragment}"

    on_disk = ROOT / resolved
    if on_disk.is_dir():
        return f"{REPO}/tree/main/{resolved}{fragment}"
    if on_disk.is_file():
        return f"{REPO}/blob/main/{resolved}{fragment}"

    FAILURES.append((source, target))
    return target


FENCE = re.compile(r"^(?P<fence>`{3,}|~{3,})(?P<info>.*)$")

# What Python-Markdown's `fenced_code` will accept as an info string. A comma is
# absent from it, which is the whole reason `normalise` exists.
LANGUAGE = re.compile(r"[\w#.+-]*")

HEADING = re.compile(r"^#{1,6}\s")


class Prepared(NamedTuple):
    """A source ready to convert, and the two counts its output must match."""

    text: str
    """The markdown, with fences normalised and links rewritten."""
    headings: int
    """`#` headings outside fences. The rendered page must have exactly this many."""
    links: int
    """`](` targets outside fences. The rendered page must have at least this many."""


def normalise(text: str, source: str, link_base: str) -> Prepared:
    """Prepare one source for conversion, line by line and fence-aware.

    Three jobs, and all three need to know where the fences are:

    **Fence info strings are cut to their first language token.** `README.md`
    opens a fence with ```` ```rust,ignore ````, which is how rustdoc spells "show
    this but do not compile it". Python-Markdown's `fenced_code` info-string
    pattern has no comma in it, so that opener does not register as a fence at
    all -- it degrades to an inline `<code>`, and its *closing* fence then opens a
    block that runs to the next fence in the file. In this corpus that silently
    ate 116 lines of `README.md`, a whole `##` section included.
    `pymdownx.superfences` rejects the same info string, so the fix belongs here
    rather than in a different extension -- and here rather than in the markdown,
    because the source is correct rustdoc that GitHub renders properly and
    `ignore` is a directive with no meaning in HTML, where only the highlight
    class matters.

    **Links are rewritten outside fences only.** A `](…)` inside a bash or Rust
    block is code; resolving it against the repository would both corrupt the
    snippet and invent a broken-link failure out of a comment.

    **Headings and links are counted**, so [`check_converted`] can hold the
    rendered page against them. That count is the only thing that would have
    caught the fence bug, which every other check in this file passed straight
    through.
    """
    out: list[str] = []
    headings = 0
    links = 0
    fence: str | None = None
    for line in text.split("\n"):
        match = FENCE.match(line)
        if fence is None and match:
            fence = match["fence"]
            language = LANGUAGE.match(match["info"].strip()).group()
            out.append(f"{fence}{language}")
            continue
        if fence is not None:
            # A closing fence is the same character, at least as long, and bare.
            closes = (
                match
                and match["fence"][0] == fence[0]
                and len(match["fence"]) >= len(fence)
                and not match["info"].strip()
            )
            if closes:
                fence = None
            out.append(line)
            continue
        headings += bool(HEADING.match(line))
        links += line.count("](")
        out.append(rewrite_links(line, source, link_base))
    if fence is not None:
        sys.exit(f"{source} has an unclosed `{fence}` fence")
    return Prepared("\n".join(out), headings, links)


def rewrite_links(text: str, source: str, link_base: str) -> str:
    """Rewrite every inline markdown link and image target in one line.

    A scanner rather than a regex, because a target can carry balanced
    parentheses -- `docs/experiments/p-60.csv` does not, but a heading anchor
    like `#the-fix-(m-348)` would, and a regex that stops at the first `)`
    truncates it into a broken link that looks deliberate.
    """
    out: list[str] = []
    at = 0
    while at < len(text):
        open_bracket = text.find("](", at)
        if open_bracket == -1:
            out.append(text[at:])
            break
        cursor = open_bracket + 2
        depth = 1
        while cursor < len(text) and depth:
            if text[cursor] == "(":
                depth += 1
            elif text[cursor] == ")":
                depth -= 1
            cursor += 1
        if depth:
            out.append(text[at:])
            break
        out.append(text[at : open_bracket + 2])
        out.append(rewrite(text[open_bracket + 2 : cursor - 1], source, link_base))
        out.append(")")
        at = cursor
    return "".join(out)


RENDERED_HEADING = re.compile(r"<h[1-6][ >]")
RENDERED_TARGET = re.compile(r"\b(?:href|src)=\"")


def check_converted(body: str, prepared: Prepared, source: str) -> None:
    """Hold the rendered page against the counts [`normalise`] took of its source.

    **A count, not a pattern**, and the ` ```rust,ignore ` bug is why. It turned
    116 lines of `README.md` into one code block, and every other check in this
    file passed straight through: no hotlink survived, every image was lazy, every
    link resolved, because none of them ask whether a page became HTML at all. The
    first version of this function searched the rendered prose for markdown that
    had survived -- and could not see this defect either, because it stripped
    `<pre>` regions first and the swallowed prose was *inside* the `<pre>`.

    Headings are an equality: markdown that got eaten takes its `##` with it, and
    a page that grew one cannot have come from this source. Link targets are a
    floor rather than an equality, because autolinks and the lazy-image rewrite
    both add `href`/`src` attributes that no `](` in the source accounts for.
    """
    headings = len(RENDERED_HEADING.findall(body))
    targets = len(RENDERED_TARGET.findall(body))
    if headings != prepared.headings:
        sys.exit(
            f"{source}: {prepared.headings} `#` headings in the source but "
            f"{headings} in the rendered page -- markdown was swallowed"
        )
    if targets < prepared.links:
        sys.exit(
            f"{source}: {prepared.links} `](` targets in the source but only "
            f"{targets} href/src in the rendered page -- markdown was swallowed"
        )


def title_of(text: str, source: str) -> str:
    """The page title: the source's first `# ` heading."""
    for line in text.splitlines():
        if line.startswith("# "):
            return line.removeprefix("# ").strip()
    sys.exit(f"{source} has no `# ` heading to take a title from")


def render(template: str, source: str, output: str, link_base: str) -> int:
    """Render one source into `web/dist`, and return the bytes written."""
    text = (ROOT / source).read_text(encoding="utf-8")
    title = title_of(text, source)
    prepared = normalise(text, source, link_base)
    body = markdown.Markdown(extensions=EXTENSIONS).convert(prepared.text)
    check_converted(body, prepared, source)
    # `DEMOS.md` carries 34 GIF references and `gameplay.md` 10, so without this
    # one substitution a single page load pulls tens of MB against a 100 GB/month
    # soft limit. The sources emit no raw HTML at all, so there is no
    # pre-existing `<img` for this to corrupt.
    body = body.replace("<img ", '<img loading="lazy" ')
    page = template.replace("{{title}}", title).replace("{{body}}", body)
    destination = OUT / output
    destination.write_text(page, encoding="utf-8")
    return len(page.encode("utf-8"))


def tree_size(path: Path) -> int:
    """Bytes under a directory, following no links."""
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    template = (ROOT / "web" / "page.html").read_text(encoding="utf-8")
    template = template.replace("{{nav}}", nav())

    total = 0
    for source, output, _, link_base in PAGES:
        written = render(template, source, output, link_base)
        total += written
        print(f"  {output:<24} {written / 1024:>8.1f} KiB  from {source}")

    if FAILURES:
        for source, target in FAILURES:
            print(f"broken link  {source}: {target}", file=sys.stderr)
        sys.exit(f"{len(FAILURES)} link(s) resolve to nothing")

    for source_dir, dest in ASSET_DIRS:
        shutil.copytree(ROOT / source_dir, OUT / dest, dirs_exist_ok=True)
        size = tree_size(OUT / dest)
        total += size
        print(f"  {dest + '/':<24} {size / 1024 / 1024:>8.1f} MiB  from {source_dir}")

    for name in ("style.css", "play.html"):
        shutil.copy2(ROOT / "web" / name, OUT / name)
        size = (OUT / name).stat().st_size
        total += size
        print(f"  {name:<24} {size / 1024:>8.1f} KiB  from web/{name}")

    print(f"  {'total':<24} {total / 1024 / 1024:>8.1f} MiB  in {OUT}")


if __name__ == "__main__":
    main()
