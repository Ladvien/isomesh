"""Trim uniform background from the bottom of a demo capture.

The examples render into a fixed window and the HUD is top-anchored, so a capture
of a small scene leaves a third of the frame empty. Crops to the content plus a
margin, top-anchored, so the HUD and the meshes are both kept.
"""
import sys
from PIL import Image

def trim(path, out, margin=72):
    im = Image.open(path).convert("RGB")
    w, h = im.size
    bg = im.getpixel((w - 4, h - 4))          # a corner that is always empty
    px = im.load()
    last = 0
    for y in range(h):
        row_has_content = False
        for x in range(0, w, 2):              # every 2nd column
            p = px[x, y]
            # A low threshold on purpose: a dark wireframe edge against a dark
            # background is exactly the content most likely to be clipped, and the
            # first version of this script clipped it.
            if abs(p[0]-bg[0]) + abs(p[1]-bg[1]) + abs(p[2]-bg[2]) > 6:
                row_has_content = True
                break
        if row_has_content:
            last = y
    bottom = min(h, last + margin)
    im.crop((0, 0, w, bottom)).save(out)
    print(f"{path}: {w}x{h} -> {w}x{bottom}  (content ends {last})")

for a in sys.argv[1:]:
    trim(a, a)
