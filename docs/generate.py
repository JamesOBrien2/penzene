"""Writes the reference pages that come from the code, so they can't drift from it:

- keys.md: Help → Keyboard Shortcuts (the table in src/MainWindow.cpp)
- python-api.md: the penzene module's docstrings

python docs/generate.py            # rewrite them (needs the built module on PYTHONPATH)
python docs/generate.py --check    # exit 1 if they're out of date (a ctest)

Both also fail if a public name's docstring doesn't end with "(since X.Y)".
"""
import html
import inspect
import pathlib
import re
import sys

DOCS = pathlib.Path(__file__).parent
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def keys_page():
    src = (DOCS.parent / "src" / "MainWindow.cpp").read_text(encoding="utf-8")
    table = re.search(r'<table cellspacing="5">(.*?)</table>', src, re.S).group(1)
    out = ["# Keyboard shortcuts", "",
           "Generated from Help → Keyboard Shortcuts in the app. Point at an atom or bond to make it the",
           "*hotspot*, then type.", ""]
    for row in re.findall(r"<tr>(.*?)</tr>", table, re.S):
        if "<th" in row:
            out += ["", "## " + text(row), "", "| Keys | |", "|---|---|"]
        else:
            cells = re.findall(r"<td>(.*?)</td>", row, re.S)
            out.append("| " + " | ".join(text(c, code=(i == 0)) for i, c in enumerate(cells)) + " |")
    return "\n".join(out) + "\n"


def text(cell, code=False):
    cell = re.sub(r"<sub>(.*?)</sub>", lambda m: m.group(1).translate(SUB), cell)
    cell = re.sub(r"<b>(.*?)</b>", (lambda m: f"`{m.group(1)}`") if code else r"**\1**", cell)
    return html.unescape(re.sub(r"<[^>]+>", "", cell)).replace("|", "\\|").strip()


def python_page():
    import penzene as pz
    out = ["# Python API reference", "",
           "Generated from the `penzene` module's docstrings. A guide with examples: [Python](python.md).", "",
           "## Functions", ""]
    for name in ("from_smiles", "read", "from_json", "drawing_styles"):
        out += member(pz, name)
    for cls in (pz.Document, pz.Atom, pz.Bond):
        out += ["", f"## {cls.__name__}", "", (cls.__doc__ or "").strip(), ""]
        for name in sorted(public(cls), key=lambda n: (n != "__init__", n)):
            out += member(cls, name, cls.__name__)
    return "\n".join(out) + "\n"


def public(cls):
    """A class's public names, with __init__ if Python code can call it (Atom and Bond have none)."""
    return [n for n, v in vars(cls).items()
            if not n.startswith("_") or (n == "__init__" and not isinstance(v, type(object.__init__)))]


def member(owner, name, prefix=None):
    obj = getattr(owner, name)
    doc = (getattr(obj, "__doc__", None) or "").strip()
    title = f"{prefix}.{name}" if prefix else name
    if isinstance(inspect.getattr_static(owner, name), property):
        return [f"### `{title}`", "", doc or "", ""]
    sig, _, rest = doc.partition("\n")  # nanobind puts the signature on the first line
    sig = sig.replace("penzene._penzene.", "")
    return [f"### `{sig if sig.startswith(name) else title}`", "", rest.strip(), ""]


def undated():
    """Public names whose docstring doesn't end with "(since X.Y)", the release that added them."""
    import penzene as pz
    names = [(pz, n) for n in pz.__all__ if n != "__version__"]
    names += [(c, n) for c in (pz.Document, pz.Atom, pz.Bond) for n in public(c)]
    return [f"{getattr(o, '__name__', 'penzene')}.{n}" for o, n in names
            if not re.search(r"\(since \d+\.\d+\)$", (getattr(o, n).__doc__ or "").strip())]


def main():
    if missing := undated():
        sys.exit(f"no \"(since X.Y)\" at the end of the docstring: {', '.join(missing)}")
    pages = {"keys.md": keys_page(), "python-api.md": python_page()}
    stale = [p for p, body in pages.items() if not (DOCS / p).exists() or (DOCS / p).read_text(encoding="utf-8") != body]
    if "--check" in sys.argv:
        if stale:
            sys.exit(f"out of date: {', '.join(stale)} (run python docs/generate.py)")
        return
    for p in stale:
        (DOCS / p).write_text(pages[p], encoding="utf-8")
        print("wrote", p)


if __name__ == "__main__":
    main()
