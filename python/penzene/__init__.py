"""Penzene: 2D chemical structures, drawn by the same engine as the app.

    import penzene as pz
    doc = pz.from_smiles("CC(=O)Oc1ccccc1C(=O)O")   # aspirin
    doc.export("aspirin.svg")
"""
import itertools as _itertools
import os as _os
import re as _re

# A wheel bundles Qt's plugins next to this file; point Qt there before it loads.
_plugins = _os.path.join(_os.path.dirname(__file__), "plugins")
if _os.path.isdir(_plugins):
    _os.environ.setdefault("QT_PLUGIN_PATH", _plugins)

from ._penzene import Atom, Bond, Document, __version__, drawing_styles, from_json, from_smiles, read  # noqa: E402

__all__ = ["Atom", "Bond", "Document", "__version__", "combinations", "drawing_styles", "from_json", "from_smiles",
           "read", "write_library"]


def combinations(sites):
    """combinations(sites: dict[str, list[str]]) -> Iterator[dict[str, str]]

    Every assignment that takes one fragment for each site, from a list per site:
    ``combinations({"R1": ["Me", "OMe"], "R2": ["Cl", "H"]})`` gives four, for write_library. (since 2.0)
    """
    keys = list(sites)
    for values in _itertools.product(*(sites[k] for k in keys)):
        yield dict(zip(keys, values))


def write_library(scaffold, assignments, path):
    """write_library(scaffold: Document, assignments: Iterable[Mapping[str, str]], path: str) -> list[tuple[int, str]]

    An SDF file of the scaffold with its R sites filled in, one record per assignment, written as each is
    read so a long table never sits in memory. An assignment maps sites ("R1", or 1) to fragments, as
    Document.substitute takes them; its other keys (an ID, a name) are metadata. Every key is kept as a
    data field of the record. Rows from ``csv.DictReader`` work as they are, or combinations() makes
    them. Returns the assignments that made no record, as (row, reason), rows counted from 1. (since 2.0)
    """
    failures = []
    with open(path, "w", encoding="utf-8", newline="\n") as out:
        for row, assignment in enumerate(assignments, 1):
            sites = {k: v for k, v in assignment.items() if isinstance(k, int) or _re.fullmatch(r"R\d+", str(k))}
            try:
                block = scaffold.substitute(sites).to_molblock()
            except ValueError as e:
                failures.append((row, str(e)))
                continue
            if not block:
                failures.append((row, "not valid chemistry"))
                continue
            out.write(block if block.endswith("\n") else block + "\n")
            for key, value in assignment.items():
                key = f"R{key}" if isinstance(key, int) else key
                out.write(f"> <{key}>\n{value}\n\n")
            out.write("$$$$\n")
    return failures
