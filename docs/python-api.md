# Python API reference

Generated from the `penzene` module's docstrings. A guide with examples: [Python](python.md).

## Functions

### `from_smiles(smiles: str) -> Document`

A structure from SMILES, laid out in 2D. (since 0.4)

### `read(path: str) -> Document`

Open a file: .penz, MOL, ChemDraw .cdxml/.cdx, .rxn, a Penzene SVG/PNG, or every record of an SDF, .smi or .inchi file laid out as a grid; records that can't be read are left out with a UserWarning naming them. (since 0.4)

### `from_json(json: str) -> Document`

A Document from .penz JSON (as to_json writes). (since 0.4)

### `drawing_styles() -> list[str]`

The drawing style names Document.style takes. (since 0.4)

### `write_library(scaffold: Document, assignments: Iterable[Mapping[str, str]], path: str) -> list[tuple[int, str]]`

An SDF file of the scaffold with its R sites filled in, one record per assignment, written as each is
read so a long table never sits in memory. An assignment maps sites ("R1", or 1) to fragments, as
Document.substitute takes them; its other keys (an ID, a name) are metadata. Every key is kept as a
data field of the record. Rows from ``csv.DictReader`` work as they are, or combinations() makes
them. Returns the assignments that made no record, as (row, reason), rows counted from 1. (since 2.0)

### `combinations(sites: dict[str, list[str]]) -> Iterator[dict[str, str]]`

Every assignment that takes one fragment for each site, from a list per site:
``combinations({"R1": ["Me", "OMe"], "R2": ["Cl", "H"]})`` gives four, for write_library. (since 2.0)


## Document

A drawing: molecules, arrows and text, as the app holds it. (since 0.4)

### `__init__(self) -> None`

An empty drawing. (since 0.4)

### `add_atom(self, label: str = 'C', x: float = 0.0, y: float = 0.0) -> int`

Add an atom: an element, an abbreviation (OMe, Boc…) or a SMILES fragment. Returns its index. (since 0.4)

### `add_bond(self, a: int, b: int, order: int = 1) -> None`

Bond atoms a and b (order 1, 2 or 3). (since 0.4)

### `Document.atoms`

The atoms, in order (their index is what add_bond and hotkeys take) (since 0.4)

### `Document.bonds`

The bonds (since 0.4)

### `clean(self, atoms: collections.abc.Sequence[int] = []) -> None`

Lay out afresh with RDKit; with atoms, only the molecules containing them. (since 0.4)

### `Document.exact_mass`

Monoisotopic mass (since 0.4)

### `export(self, path: str) -> None`

Write .svg, .png or .pdf, exactly as the app exports. (since 0.4)

### `Document.formula`

Molecular formula, Hill order, all molecules together (since 0.4)

### `hotkeys(self, atom: int, keys: str, bond: int | None = None) -> tuple[int, int]`

Type hotkeys with this atom as the hotspot (or atom=-1, bond=i for a bond). Returns the final (atom, bond) hotspot. (since 0.4)

### `Document.mw`

Molecular weight (since 0.4)

### `Document.r_sites`

The numbers of the R sites (atoms labelled R1, R2…, or [1*] and [*:1] in SMILES), sorted. (since 2.0)

### `remove_atoms(self, atoms: collections.abc.Sequence[int]) -> None`

Remove atoms and their bonds. The atoms after each one removed move down to fill its place, so indices taken before are stale: remove everything in one call. (since 2.0)

### `remove_bonds(self, bonds: collections.abc.Sequence[int]) -> None`

Remove bonds; their atoms stay. Later bonds move down as in remove_atoms, and atom indices don't change. (since 2.0)

### `save(self, path: str) -> None`

Write .penz (full fidelity), .mol, or ChemDraw .cdxml or .cdx; images go through export(). (since 0.4)

### `set_bond_order(self, bond: int, order: int) -> None`

Make a bond single, double or triple (a wedge or hash goes with a single bond). (since 2.0)

### `set_label(self, atom: int, label: str) -> None`

Relabel an atom, as add_atom takes a label: an element ("C" undoes a label), an abbreviation or a SMILES fragment (drawn out from it, its new atoms added at the end). (since 2.0)

### `Document.style`

Drawing style: 'ACS 1996' (default), 'JDP' or 'RSC' (since 0.4)

### `substitute(self, fragments: dict) -> Document`

A copy with each R site replaced, from a dict keyed by site (1 or "R1"): an abbreviation (OMe, Ph) or SMILES whose first atom, or the atom bonded to its lone * ("*OC"), takes the site's place; "H" removes it. The scaffold keeps its drawing and wedges. ValueError names a missing or unknown site, an R site with more than one bond, or a fragment that doesn't parse. (since 2.0)

### `to_inchi(self) -> str`

Standard InChI. (since 0.4)

### `to_inchikey(self) -> str`

Standard InChIKey. (since 0.4)

### `to_json(self) -> str`

The .penz document (JSON). (since 0.4)

### `to_molblock(self) -> str`

MDL MOL, with the drawing's wedges; V3000 when coordinate bonds require it. (since 0.4)

### `to_png(self, dpi: float = 300) -> bytes`

PNG bytes at this resolution. (since 0.4)

### `to_smiles(self) -> str`

Canonical SMILES ("" if the drawing isn't valid chemistry), as CXSMILES carrying only the stereo groups (&1, or1) when there are any. (since 0.4)

### `to_svg(self) -> str`

SVG, as the app exports it (the drawing embedded, so it reopens editable). (since 0.4)


## Atom

An atom of a Document (read-only; edit through Document: set_label, remove_atoms). (since 0.4)

### `Atom.charge`

Formal charge (since 0.4)

### `Atom.isotope`

Mass number (13 for carbon-13, 2 for deuterium); 0 for natural abundance (since 1.3)

### `Atom.symbol`

Element symbol, or the abbreviation (Boc, OMe…) if it has one (since 0.4)

### `Atom.x`

Points, x right (since 0.4)

### `Atom.y`

Points, y down (since 0.4)


## Bond

A bond of a Document (read-only). (since 0.4)

### `Bond.a`

Index of the first atom (a wedge starts here) (since 0.4)

### `Bond.b`

Index of the second atom (since 0.4)

### `Bond.order`

1, 2 or 3 (since 0.4)

