"""Builds docs/_static/tutorials/caffeine-keys.gif and its last frame from the caffeine test's frames.

    PENZENE_TUTORIAL_FRAMES=/tmp/frames QT_QPA_PLATFORM=offscreen build/bin/penzene_tests "caffeine drawn*"
    pixi exec --spec pillow --spec python python docs/tutorials/make-gif.py /tmp/frames
"""
import os
import sys

from PIL import Image, ImageDraw, ImageFont

frames_dir = sys.argv[1]
out = os.path.join(os.path.dirname(__file__), "..", "_static", "tutorials")
os.makedirs(out, exist_ok=True)
with open(os.path.join(frames_dir, "frames.txt"), encoding="utf-8") as f:
    steps = [line.rstrip("\n").split("\t") for line in f if line.strip()]

font = None
for path in ("/System/Library/Fonts/Helvetica.ttc", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
    if os.path.exists(path):
        font = ImageFont.truetype(path, 15)
        break
font = font or ImageFont.load_default()

frames = []
for name, caption in steps:
    shot = Image.open(os.path.join(frames_dir, name)).convert("RGB")
    frame = Image.new("RGB", (shot.width, shot.height + 34), "white")
    frame.paste(shot, (0, 0))
    draw = ImageDraw.Draw(frame)
    draw.rectangle([0, shot.height, shot.width, shot.height + 34], fill=(238, 244, 242))
    draw.text((10, shot.height + 8), caption, fill=(20, 40, 36), font=font)
    frames.append(frame)
# Each step shows the hotspot, then the result; the first and last frames stay longer.
durations = [2500] + [700] * (len(frames) - 2) + [3500]
frames[0].save(os.path.join(out, "caffeine-keys.gif"), save_all=True, append_images=frames[1:],
               duration=durations, loop=0, optimize=True)
frames[-1].save(os.path.join(out, "caffeine-keys.png"), optimize=True)
print("wrote", os.path.join(out, "caffeine-keys.gif"))
