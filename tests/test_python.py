"""The Python module: the #47 snippet as written, plus the rest of the API."""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import struct
import zlib
import xml.etree.ElementTree as ET

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

# Atom and Bond snapshots stay read-only, with a useful route to editing (#716).
before = doc.to_json()
for obj, fields in (
    (doc.atoms[0], dict(symbol="set_label", x="set_position", y="set_position", charge="set_charge",
                        isotope="Document", label="set_label", radicals="set_radicals", partial="set_partial")),
    (doc.bonds[0], dict(a="add_bond", b="add_bond", order="set_bond_order", style="set_bond_style")),
):
    for field, hint in fields.items():
        original = getattr(obj, field)
        for value in (original, None):
            try:
                setattr(obj, field, value)
            except AttributeError as error:
                assert f"{type(obj).__name__}.{field} is read-only" in str(error), str(error)
                assert hint in str(error), str(error)
            else:
                raise AssertionError(f"{field} should be read-only")
        assert getattr(obj, field) == original
assert doc.to_json() == before

# Export options use the app's frame, while bounds exclude the requested margin (#719).
def png_info(data):
    width, height, depth, color = struct.unpack(">IIBB", data[16:26])
    assert depth == 8 and color in (2, 6)
    chunks = []
    offset = 8
    while offset < len(data):
        size = struct.unpack(">I", data[offset:offset + 4])[0]
        if data[offset + 4:offset + 8] == b"IDAT":
            chunks.append(data[offset + 8:offset + 8 + size])
        offset += size + 12
    # Every PNG filter leaves the first pixel unchanged on the first row.
    pixel = zlib.decompress(b"".join(chunks))[1:1 + (4 if color == 6 else 3)]
    return (width, height), tuple(pixel)

x0, y0, x1, y1 = doc.bounds()
assert x0 < x1 and y0 < y1
frame = [float(x) for x in ET.fromstring(doc.to_svg(margin=3, scale=2, background="white")).attrib["viewBox"].split()]
assert abs(frame[2] - (2 * (x1 - x0) + 6)) < 1e-5
assert abs(frame[3] - (2 * (y1 - y0) + 6)) < 1e-5
png = doc.to_png(72, margin=3, background="#ffffff", scale=2)
size, pixel = png_info(png)
assert size == (round(frame[2]), round(frame[3])) and pixel in ((255, 255, 255), (255, 255, 255, 255))
assert png_info(doc.to_png(72, margin=3))[1][-1] == 0
assert doc.to_svg() == doc.to_svg(margin=0, background=None, scale=1)
for ext in ("png", "svg", "pdf"):
    path = os.path.join(out, "options." + ext)
    doc.export(path, margin=3, background="white", scale=2, dpi=72)
    with open(path, "rb") as f:
        content = f.read()
    if ext == "png":
        assert content == png
    elif ext == "svg":
        assert content.decode() == doc.to_svg(margin=3, background="white", scale=2)
    else:
        assert content.startswith(b"%PDF")
mark = pz.Document()
mark.add_atom("R", generic=True), mark.set_charge(0, -1), mark.set_radicals(0, 1)
mark_bounds = mark.bounds()
mark.add_atom("O", x=200), mark.add_text(400, 0, "free text")
assert mark.bounds([0, 0]) == mark_bounds and mark.bounds()[2] > 400
assert pz.Document().bounds() == (0, 0, 0, 0)
styled_bounds = pz.from_smiles("CO")
styled_bounds.style = "RSC"
bounds = styled_bounds.bounds()
frame = [float(x) for x in ET.fromstring(styled_bounds.to_svg()).attrib["viewBox"].split()]
assert abs(frame[2] - (bounds[2] - bounds[0])) < 1e-4  # style's physical points, not model coordinates
for bad in (lambda: doc.bounds([0, 999]), lambda: doc.bounds([-1])):
    try:
        bad()
    except IndexError:
        pass
    else:
        raise AssertionError("bounds must validate atom indices")
for options in ({"margin": -1}, {"margin": float("nan")}, {"scale": 0}, {"scale": float("inf")},
                {"scale": 1e300}, {"background": "not-a-colour"}):
    for render in (doc.to_svg, doc.to_png, lambda **opts: doc.export(os.path.join(out, "invalid.png"), **opts)):
        try:
            render(**options)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid export options accepted: {options}")
assert not os.path.exists(os.path.join(out, "invalid.png"))
for dpi in (0, -1, float("inf"), 1e100):
    try:
        doc.to_png(dpi)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid dpi accepted")

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
for journal in ("Nature", "Angewandte", "SYNLETT/SYNTHESIS"):
    assert journal in pz.drawing_styles()
    doc.style = journal
    assert doc.style == journal
    assert "<svg" in doc.to_svg()
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
import warnings
with warnings.catch_warnings(record=True) as caught:  # a partial library warns, naming the record (#584)
    warnings.simplefilter("always")
    partial = pz.read(os.path.join(os.environ["PENZENE_TEST_DATA"], "broken.smi"))
assert partial.formula == "C2H6O"
assert any(issubclass(w.category, UserWarning) and "broken" in str(w.message) for w in caught), caught
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
styled = pz.from_smiles("C=C")
styled.set_bond_style(0, "double-dashed")  # keeps the double bond's order, both lines dashed
assert styled.bonds[0].order == 2 and "stroke-dasharray" in styled.to_svg()
quad = pz.from_smiles("C=C")
quad.set_bond_order(0, 4)
assert quad.bonds[0].order == 4
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

# Atom-pair bond lookup always uses current indices, including after removals (#720).
lookup = pz.from_smiles("CCCO")
assert lookup.bond_between(1, 2) == lookup.bond_between(2, 1) == 1
assert lookup.bond_between(0, 3) is None and lookup.bond_between(0, 0) is None
lookup.set_bond_order(lookup.bond_between(1, 2), 2)
assert lookup.bonds[1].order == 2
lookup.remove_atoms([0])
assert lookup.bond_between(0, 1) == 0 and lookup.bonds[0].order == 2
lookup.remove_bonds([0])
assert lookup.bond_between(0, 1) is None and lookup.bond_between(1, 2) == 0
for a, b in ((-1, 0), (0, -1), (3, 0), (0, 3)):
    try:
        lookup.bond_between(a, b)
    except IndexError:
        pass
    else:
        raise AssertionError("bond_between must validate both atom indices")

# Reaction schemes from Python use the app's arrow kinds and editable text (#718).
scheme = pz.Document()
for i, kind in enumerate(("reaction", "resonance", "equilibrium", "retro")):
    assert scheme.add_arrow(0, i * 50, 60, i * 50, kind=kind, text_above="H2O", text_below="heat") == i
scheme_page = json.loads(scheme.to_json())
scheme_page = (scheme_page.get("pages") or [scheme_page])[0]
assert [a["kind"] for a in scheme_page["arrows"]] == ["reaction", "resonance", "equilibrium", "retro"]
assert len(scheme_page["texts"]) == 8
for i in range(4):
    above, below = scheme_page["texts"][2 * i:2 * i + 2]
    assert above["text"] == "H2O" and above["y"] < i * 50
    assert below["text"] == "heat" and below["y"] > i * 50
assert pz.from_json(scheme.to_json()).to_svg() == scheme.to_svg()
assert "<path" in scheme.to_svg() and scheme.to_png().startswith(b"\x89PNG")
assert scheme.add_arrow(100, 0, 0, 50) == 4  # leftward, diagonal, and default kind
page = json.loads(scheme.to_json())
page = (page.get("pages") or [page])[0]
assert page["arrows"][-1]["kind"] == "reaction" and len(page["texts"]) == 8
for bad in (lambda: scheme.add_arrow(0, 0, 0, 0),
            lambda: scheme.add_arrow(0, 0, 20, 0, kind="unknown"),
            lambda: scheme.add_arrow(float("nan"), 0, 20, 0),
            lambda: scheme.add_arrow(0, 0, float("inf"), 0)):
    before = scheme.to_json()
    try:
        bad()
    except ValueError:
        assert scheme.to_json() == before
    else:
        raise AssertionError("invalid arrow must fail without changing the drawing")

# Atom transforms (#717): selected atoms move once, with bonds and stereo intact.
moved = pz.Document()
for x, y in ((0, 0), (10, 0), (20, 10)):
    moved.add_atom("C", x, y)
moved.add_bond(0, 1, 2), moved.add_bond(1, 2)
moved.move_atoms([0, 1, 0], 2, 3)
assert [(a.x, a.y) for a in moved.atoms] == [(2, 3), (12, 3), (20, 10)]
moved.set_position(2, 22, 3)
moved.rotate(degrees=90, about=(2, 3))
for a, expected in zip(moved.atoms, ((2, 3), (2, 13), (2, 23))):
    assert abs(a.x - expected[0]) < 1e-8 and abs(a.y - expected[1]) < 1e-8
assert [(b.a, b.b, b.order) for b in moved.bonds] == [(0, 1, 2), (1, 2, 1)]
moved.rotate([0, 1], degrees=180)  # default pivot is this selection's mean
assert abs(moved.atoms[0].y - 13) < 1e-8 and abs(moved.atoms[1].y - 3) < 1e-8
assert abs(moved.atoms[2].y - 23) < 1e-8
mirror_data = json.loads(moved.to_json())
(mirror_data.get("pages") or [mirror_data])[0]["bonds"][0]["position"] = "left"
mirrored = pz.from_json(json.dumps(mirror_data))
mirrored.flip([0, 0, 1])
mirror_page = json.loads(mirrored.to_json())
assert (mirror_page.get("pages") or [mirror_page])[0]["bonds"][0]["position"] == "right"
mirrored.flip([0])  # moving just one end does not swap the explicit side
mirror_page = json.loads(mirrored.to_json())
assert (mirror_page.get("pages") or [mirror_page])[0]["bonds"][0]["position"] == "right"
chiral_move = pz.from_smiles("C[C@H](N)C(=O)O")
original_smiles = chiral_move.to_smiles()
chiral_move.rotate(degrees=37)
assert chiral_move.to_smiles() == original_smiles
chiral_move.flip()
assert chiral_move.to_smiles() == pz.from_smiles("C[C@@H](N)C(=O)O").to_smiles()
chiral_move.flip(horizontal=False)
assert chiral_move.to_smiles() == original_smiles
for bad, kind in ((lambda: moved.move_atoms([0, 99], 1, 2), IndexError),
                  (lambda: moved.set_position(-1, 0, 0), IndexError),
                  (lambda: moved.flip([99]), IndexError),
                  (lambda: moved.rotate([0, 99], degrees=90), IndexError),
                  (lambda: moved.move_atoms([], float("nan"), 0), ValueError),
                  (lambda: moved.set_position(0, float("inf"), 0), ValueError),
                  (lambda: moved.rotate(degrees=float("inf")), ValueError),
                  (lambda: moved.rotate(degrees=1, about=(0, float("nan"))), ValueError)):
    before = moved.to_json()
    try:
        bad()
    except kind:
        assert moved.to_json() == before
    else:
        raise AssertionError("invalid transform must fail without changing the drawing")
empty = pz.Document()
empty.move_atoms([], 1, 2), empty.flip(), empty.rotate(degrees=90)

# Mechanism marks from Python (#659): a three-centre TS for R-NC to R-CN, with no JSON editing.
ts = pz.Document()
r, n, c = ts.add_atom("R", generic=True), ts.add_atom("N", 16, -7.6), ts.add_atom("C", 16, 7.6)
assert ts.atoms[r].label == "R"
ts.add_bond(r, n), ts.add_bond(n, c, 3), ts.add_bond(r, c)
ts.set_bond_style(0, "partial"), ts.set_bond_style(2, "partial")
ts.set_partial(r, 1), ts.set_partial(c, -1)
ts.show_symbol(c)
assert ts.add_bracket([c, n, r, r], superscript="\u2021") == 0
assert ts.add_text(0, 30, "TS1") == 0
assert [b.style for b in ts.bonds] == ["partial", "none", "partial"]
assert (ts.atoms[r].partial, ts.atoms[c].partial, ts.atoms[c].label, ts.atoms[n].label) == (1, -1, "C", "")
saved = json.loads(ts.to_json())
page = (saved.get("pages") or [saved])[0]
assert page["brackets"] == [{"atoms": [0, 1, 2], "superscript": "\u2021"}] and page["texts"][0]["text"] == "TS1"
ts.show_symbol(c, False)
assert ts.atoms[c].label == ""
ion = pz.from_smiles("C[CH2]")
ion.set_charge(1, -1), ion.set_radicals(0, 1)
assert (ion.atoms[1].charge, ion.atoms[0].radicals) == (-1, 1)
for bad, kind in ((lambda: ts.set_bond_style(1, "wedge"), ValueError), (lambda: ts.set_bond_style(0, "squiggly"), ValueError),
                  (lambda: ts.set_charge(9, 1), IndexError), (lambda: ts.set_radicals(0, 3), ValueError),
                  (lambda: ts.show_symbol(r), ValueError), (lambda: ts.add_bracket([]), ValueError)):
    try:
        bad()
    except kind:
        pass
    else:
        raise AssertionError("should raise " + kind.__name__)

# R-group libraries (#432): numbered wildcards are R sites, filled from fragments into an SDF.
scaffold = pz.from_smiles("[1*]c1ccc([*:2])cc1")
assert scaffold.r_sites == [1, 2]
filled = scaffold.substitute({"R1": "OC", 2: "*C(F)(F)F"})
assert filled.to_smiles() == pz.from_smiles("COc1ccc(C(F)(F)F)cc1").to_smiles()
assert scaffold.substitute({1: "Cl", 2: "H"}).to_smiles() == pz.from_smiles("Clc1ccccc1").to_smiles()
assert scaffold.substitute({1: "OMe", 2: "Ph"}).formula == "C13H12O"
chiral = pz.from_smiles("C[C@H]([1*])N")  # the scaffold's stereocentre survives the substitution
assert chiral.substitute({1: "O"}).to_smiles() == pz.from_smiles("C[C@H](O)N").to_smiles()
for bad, why in [({1: "OC"}, "R2"), ({1: "OC", 2: "C", 3: "C"}, "R3"), ({1: "OC", 2: "C1CC"}, "R2"),
                 ({1: "OC", 2: "*C*"}, "R2"), ({"X": "C"}, "X")]:
    try:
        scaffold.substitute(bad)
        raise AssertionError(f"{bad} must raise")
    except ValueError as e:
        assert why in str(e), (bad, e)
rows = list(pz.combinations({"R1": ["Cl", "OC"], "R2": ["H", "nonsense("]}))
assert len(rows) == 4 and rows[0] == {"R1": "Cl", "R2": "H"}
rows = [dict(r, ID=f"cpd-{k}") for k, r in enumerate(rows, 1)]
library = os.path.join(out, "library.sdf")
assert pz.write_library(scaffold, iter(rows), library) == [(2, "R2: not a fragment: nonsense("), (4, "R2: not a fragment: nonsense(")]
records = pz.read(library)  # the two good rows, laid out as a grid
assert records.to_smiles().count(".") == 1
with open(library) as f:
    sdf = f.read()
assert sdf.count("$$$$") == 2 and "> <ID>\ncpd-3\n" in sdf and "> <R1>\nOC\n" in sdf

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
    with open(os.path.join(root, "docs", "_static", "theme.schema.json")) as f:
        theme_schema = json.load(f)
    jsonschema.Draft202012Validator.check_schema(theme_schema)
    with open(os.path.join(root, "docs", "_static", "theme-example.json")) as f:
        jsonschema.validate(json.load(f), theme_schema)
    examples = sorted(glob.glob(os.path.join(os.environ["PENZENE_TEST_DATA"], "penz", "v*.penz")))
    assert len(examples) >= 8
    for path in examples:
        with open(path) as f:
            jsonschema.validate(json.load(f), schema)
        jsonschema.validate(json.loads(pz.read(path).to_json()), schema)
