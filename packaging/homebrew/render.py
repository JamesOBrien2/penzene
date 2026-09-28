#!/usr/bin/env python3
"""Fill penzene.rb.in for a release: render.py 1.4.0 > penzene.rb"""
import hashlib, json, pathlib, sys, urllib.request

version = sys.argv[1].removeprefix("v")
api = f"https://api.github.com/repos/JamesOBrien2/penzene/releases/tags/v{version}"
assets = {a["name"]: a for a in json.load(urllib.request.urlopen(api))["assets"]}


def sha256(name):
    digest = assets[name].get("digest") or ""  # "sha256:...", when GitHub has one
    if digest.startswith("sha256:"):
        return digest.removeprefix("sha256:")
    with urllib.request.urlopen(assets[name]["browser_download_url"]) as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


text = (pathlib.Path(__file__).with_name("penzene.rb.in").read_text()
        .replace("@VERSION@", version)
        .replace("@SHA256_ARM64@", sha256("penzene-macos-arm64.dmg"))
        .replace("@SHA256_X86_64@", sha256("penzene-macos-x86_64.dmg")))
sys.stdout.write(text)
