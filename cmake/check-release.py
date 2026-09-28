"""Checks a published release: python3 cmake/check-release.py v1.4.0

Its GitHub assets, its wheels on PyPI, and that Read the Docs "stable" was built from its tag.
Prints one line per gap and exits 1 if there are any.
"""
import json
import sys
import urllib.error
import urllib.request

REPO = "https://api.github.com/repos/JamesOBrien2/penzene"
APPS = {"penzene-macos-arm64.dmg", "penzene-macos-x86_64.dmg", "penzene-windows-x64.zip",
        "penzene-windows-x64-setup.exe", "penzene-linux-x86_64.AppImage"}
WHEELS = 4


def get(url):
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        raise LookupError(f"{url}: HTTP {e.code}") from None


tag = sys.argv[1]
version = tag.removeprefix("v")
gaps = []

try:
    names = {a["name"] for a in get(f"{REPO}/releases/tags/{tag}")["assets"]}
    wheels = {n for n in names if n.startswith(f"penzene-{version}-") and n.endswith(".whl")}
    gaps += [f"GitHub release: missing {n}" for n in sorted(APPS - names)]
    if len(wheels) != WHEELS:
        gaps.append(f"GitHub release: {len(wheels)} wheels, expected {WHEELS}")
    gaps += [f"GitHub release: unexpected asset {n}" for n in sorted(names - APPS - wheels)]
except LookupError as e:
    gaps.append(f"GitHub release: {e}")

try:
    files = get(f"https://pypi.org/pypi/penzene/{version}/json")["urls"]
    n = sum(f["filename"].endswith(".whl") for f in files)
    if n != WHEELS:
        gaps.append(f"PyPI: {n} wheels for {version}, expected {WHEELS}")
except LookupError as e:
    gaps.append(f"PyPI: {e}")

try:
    obj = get(f"{REPO}/git/ref/tags/{tag}")["object"]
    if obj["type"] == "tag":  # annotated: the commit is one step further
        obj = get(obj["url"])["object"]
    built = get("https://readthedocs.org/api/v3/projects/penzene/versions/stable/")["identifier"]
    if built != obj["sha"]:
        gaps.append(f"Read the Docs: stable is built from {built[:12]}, {tag} is {obj['sha'][:12]}")
except LookupError as e:
    gaps.append(f"Read the Docs: {e}")

for g in gaps:
    print(g)
if gaps:
    sys.exit(1)
print(f"{tag}: {len(APPS) + WHEELS} assets, {WHEELS} wheels on PyPI, Read the Docs stable at the tag")
