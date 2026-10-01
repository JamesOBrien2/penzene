"""The Python module: the #47 snippet as written, plus the rest of the API."""
import importlib.util
import os
import subprocess
import sys
import tempfile

import penzene as pz

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
out = tempfile.mkdtemp()

doc = pz.from_smiles(ASPIRIN)
doc.clean()
doc.save(os.path.join(out, "aspirin.penz"))
doc.export(os.path.join(out, "aspirin.svg"))
assert "<svg" in doc._repr_svg_()
assert "<path" in doc.to_svg(), "labels must render as glyph outlines"
assert doc.to_png(dpi=72)[:8] == b"\x89PNG\r\n\x1a\n"
assert doc.formula == "C9H8O4" and abs(doc.mw - 180.16) < 0.01
assert doc.to_inchikey() == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
assert len(doc.atoms) == 13 and doc.atoms[0].symbol == "C"

again = pz.read(os.path.join(out, "aspirin.penz"))
assert again.to_smiles() == doc.to_smiles()
assert pz.from_json(doc.to_json()).to_smiles() == doc.to_smiles()
for ext in ("cdxml", "cdx"):  # ChemDraw, both ways
    doc.save(os.path.join(out, "aspirin." + ext))
    assert pz.read(os.path.join(out, "aspirin." + ext)).to_smiles() == doc.to_smiles(), ext
assert pz.read(os.path.join(os.environ["PENZENE_TEST_DATA"], "aspirin.mol")).formula == "C9H8O4"
try:  # save never writes MOL text into an image (#327)
    doc.save(os.path.join(out, "aspirin.png"))
except ValueError:
    pass
else:
    raise AssertionError("save to .png should raise")
assert not os.path.exists(os.path.join(out, "aspirin.png"))

# The hotkey builder: ChemDraw's cheat-sheet dipeptide, 42n152o from H2N-CH3.
m = pz.Document()
a = m.add_atom("N")
b = m.add_atom("C", x=14.4)
m.add_bond(a, b)
m.hotkeys(b, "42n152o")
assert m.formula == "C6H12N2O3", m.formula  # Ala-Ala

m2 = pz.Document()
a = m2.add_atom("N")
m2.hotkeys(a, "42n152o")  # the issue's snippet, from a lone N: must just work
assert m2.to_smiles()

boc = pz.Document()
boc.add_atom("Boc")
assert boc.atoms[0].symbol == "Boc" and boc.formula == "C5H10O2"

doc.style = "RSC"
assert doc.style == "RSC" and "RSC" in pz.drawing_styles()
for bad in (lambda: pz.from_smiles("C1CC"), lambda: m.add_atom("notachem!!"), lambda: setattr(doc, "style", "Comic")):
    try:
        bad()
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")
try:
    pz.read(os.path.join(out, "missing.penz"))
except FileNotFoundError as e:
    assert isinstance(e, ValueError), "a missing file is still a ValueError"
else:
    raise AssertionError("expected FileNotFoundError")
try:
    pz.Document().style = "nope"
except ValueError as e:
    assert "'RSC'" in str(e), e
else:
    raise AssertionError("expected ValueError")
# Editing an imported drawing in place (#430): relabel, change a bond, remove atoms and bonds.
ed = pz.from_smiles("CC(=O)OC")  # methyl acetate
ed.set_label(4, "Et")
assert ed.atoms[4].symbol == "Et" and ed.formula == "C4H8O2"
ed.set_label(4, "C")  # back to a plain carbon
assert ed.atoms[4].symbol == "C" and ed.formula == "C3H6O2"
double = next(i for i, b in enumerate(ed.bonds) if b.order == 2)
ed.set_bond_order(double, 1)
assert ed.bonds[double].order == 1 and ed.formula == "C3H8O2"
gone = pz.from_smiles("CC(=O)OC")
gone.remove_atoms([3, 1, 3])  # indices from before the call; duplicates are fine
assert [a.symbol for a in gone.atoms] == ["C", "O", "C"] and len(gone.bonds) == 0
ed.remove_atoms([4])  # the methyl
assert len(ed.atoms) == 4 and all(b.a < 4 and b.b < 4 for b in ed.bonds)
pairs = [(b.a, b.b) for b in ed.bonds]
ed.remove_bonds([1, 0])
assert [(b.a, b.b) for b in ed.bonds] == pairs[2:] and len(ed.atoms) == 4  # those two go; the atoms stay
for bad in (lambda: ed.set_label(99, "N"), lambda: ed.remove_atoms([99]), lambda: ed.set_bond_order(99, 1)):
    try:
        bad()
        raise AssertionError("an index out of range must raise")
    except IndexError:
        pass
try:
    ed.set_label(0, "notachem!!")
    raise AssertionError("an unknown label must raise")
except ValueError:
    assert ed.atoms[0].symbol == "C"  # and leave the atom as it was

print("python ok", pz.__version__)

# Type stubs: in a development build, the committed ones match the module (regenerate with
# `pixi run stubs`); from an installed wheel, they ship with it.
if importlib.util.find_spec("nanobind"):
    stub = os.path.join(out, "_penzene.pyi")
    subprocess.run([sys.executable, "-m", "nanobind.stubgen", "-q", "-m", "penzene._penzene", "-o", stub], check=True)
    package = os.path.join(os.path.dirname(os.environ["PENZENE_TEST_DATA"]), "..", "python", "penzene")
    with open(stub) as fresh, open(os.path.join(package, "_penzene.pyi")) as committed:
        assert fresh.read() == committed.read(), "type stubs are out of date: pixi run stubs"
else:
    package = os.path.dirname(pz.__file__)
    assert os.path.exists(os.path.join(package, "_penzene.pyi")), "the wheel has no type stubs"
assert os.path.exists(os.path.join(package, "py.typed"))
assert set(pz.__all__) == {n for n in dir(pz) if not n.startswith("_")} | {"__version__"}, "__all__ is the public API"

# .penz: every earlier release's example file, and what this release writes from each, fit the
# published schema (docs/_static/penz.schema.json). Skipped where jsonschema isn't installed.
if importlib.util.find_spec("jsonschema"):
    import glob
    import json

    import jsonschema

    root = os.path.join(os.path.dirname(os.environ["PENZENE_TEST_DATA"]), "..")
    with open(os.path.join(root, "docs", "_static", "penz.schema.json")) as f:
        schema = json.load(f)
    jsonschema.Draft202012Validator.check_schema(schema)
    examples = sorted(glob.glob(os.path.join(os.environ["PENZENE_TEST_DATA"], "penz", "v*.penz")))
    assert len(examples) >= 8
    for path in examples:
        with open(path) as f:
            jsonschema.validate(json.load(f), schema)
        jsonschema.validate(json.loads(pz.read(path).to_json()), schema)
