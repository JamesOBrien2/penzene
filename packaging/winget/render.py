#!/usr/bin/env python3
"""Fill the winget manifests for a release: render.py 1.4.0 [outdir]

Writes manifests/j/JamesOBrien2/Penzene/<version>/ under outdir (default: the current
directory), the layout microsoft/winget-pkgs uses.
"""
import hashlib, json, pathlib, sys, urllib.request

version = sys.argv[1].removeprefix("v")
out = pathlib.Path(sys.argv[2] if len(sys.argv) > 2 else ".") / "manifests/j/JamesOBrien2/Penzene" / version
api = f"https://api.github.com/repos/JamesOBrien2/penzene/releases/tags/v{version}"
release = json.load(urllib.request.urlopen(api))
asset = next(a for a in release["assets"] if a["name"] == "penzene-windows-x64-setup.exe")

digest = asset.get("digest") or ""  # "sha256:...", when GitHub has one
if digest.startswith("sha256:"):
    sha = digest.removeprefix("sha256:")
else:
    with urllib.request.urlopen(asset["browser_download_url"]) as f:
        sha = hashlib.file_digest(f, "sha256").hexdigest()

out.mkdir(parents=True, exist_ok=True)
for template in pathlib.Path(__file__).parent.glob("*.yaml.in"):
    text = (template.read_text().replace("@VERSION@", version).replace("@SHA256@", sha.upper())
            .replace("@DATE@", release["published_at"][:10]))
    (out / template.name.removesuffix(".in")).write_text(text)
print(out)
