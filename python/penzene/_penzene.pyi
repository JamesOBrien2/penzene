"""Penzene: 2D chemical structures, drawn by the same engine as the app."""

from collections.abc import Sequence


class Atom:
    """
    An atom of a Document (read-only; edit through Document: set_label, remove_atoms). (since 0.4)
    """

    @property
    def symbol(self) -> str:
        """
        Element symbol, or the abbreviation (Boc, OMe…) if it has one (since 0.4)
        """

    @property
    def x(self) -> float:
        """Points, x right (since 0.4)"""

    @property
    def y(self) -> float:
        """Points, y down (since 0.4)"""

    @property
    def charge(self) -> int:
        """Formal charge (since 0.4)"""

    @property
    def isotope(self) -> int:
        """
        Mass number (13 for carbon-13, 2 for deuterium); 0 for natural abundance (since 1.3)
        """

    def __repr__(self) -> str: ...

class Bond:
    """A bond of a Document (read-only). (since 0.4)"""

    @property
    def a(self) -> int:
        """Index of the first atom (a wedge starts here) (since 0.4)"""

    @property
    def b(self) -> int:
        """Index of the second atom (since 0.4)"""

    @property
    def order(self) -> int:
        """1, 2 or 3 (since 0.4)"""

    def __repr__(self) -> str: ...

class Document:
    """
    A drawing: molecules, arrows and text, as the app holds it. (since 0.4)
    """

    def __init__(self) -> None:
        """An empty drawing. (since 0.4)"""

    @property
    def atoms(self) -> list[Atom]:
        """
        The atoms, in order (their index is what add_bond and hotkeys take) (since 0.4)
        """

    @property
    def bonds(self) -> list[Bond]:
        """The bonds (since 0.4)"""

    @property
    def style(self) -> str:
        """Drawing style: 'ACS 1996' (default), 'JDP' or 'RSC' (since 0.4)"""

    @style.setter
    def style(self, arg: str, /) -> None: ...

    def add_atom(self, label: str = 'C', x: float = 0.0, y: float = 0.0) -> int:
        """
        Add an atom: an element, an abbreviation (OMe, Boc…) or a SMILES fragment. Returns its index. (since 0.4)
        """

    def add_bond(self, a: int, b: int, order: int = 1) -> None:
        """Bond atoms a and b (order 1, 2 or 3). (since 0.4)"""

    def set_label(self, atom: int, label: str) -> None:
        """
        Relabel an atom, as add_atom takes a label: an element ("C" undoes a label), an abbreviation or a SMILES fragment (drawn out from it, its new atoms added at the end). (since 2.0)
        """

    def set_bond_order(self, bond: int, order: int) -> None:
        """
        Make a bond single, double or triple (a wedge or hash goes with a single bond). (since 2.0)
        """

    def remove_atoms(self, atoms: Sequence[int]) -> None:
        """
        Remove atoms and their bonds. The atoms after each one removed move down to fill its place, so indices taken before are stale: remove everything in one call. (since 2.0)
        """

    def remove_bonds(self, bonds: Sequence[int]) -> None:
        """
        Remove bonds; their atoms stay. Later bonds move down as in remove_atoms, and atom indices don't change. (since 2.0)
        """

    @property
    def r_sites(self) -> list[int]:
        """
        The numbers of the R sites (atoms labelled R1, R2…, or [1*] and [*:1] in SMILES), sorted. (since 2.0)
        """

    def substitute(self, fragments: dict) -> Document:
        """
        A copy with each R site replaced, from a dict keyed by site (1 or "R1"): an abbreviation (OMe, Ph) or SMILES whose first atom, or the atom bonded to its lone * ("*OC"), takes the site's place; "H" removes it. The scaffold keeps its drawing and wedges. ValueError names a missing or unknown site, an R site with more than one bond, or a fragment that doesn't parse. (since 2.0)
        """

    def hotkeys(self, atom: int, keys: str, bond: int | None = None) -> tuple[int, int]:
        """
        Type hotkeys with this atom as the hotspot (or atom=-1, bond=i for a bond). Returns the final (atom, bond) hotspot. (since 0.4)
        """

    def clean(self, atoms: Sequence[int] = []) -> None:
        """
        Lay out afresh with RDKit; with atoms, only the molecules containing them. (since 0.4)
        """

    def to_smiles(self) -> str:
        """
        Canonical SMILES ("" if the drawing isn't valid chemistry), as CXSMILES carrying only the stereo groups (&1, or1) when there are any. (since 0.4)
        """

    def to_molblock(self) -> str:
        """
        MDL MOL, with the drawing's wedges; V3000 when coordinate bonds require it. (since 0.4)
        """

    def to_inchi(self) -> str:
        """Standard InChI. (since 0.4)"""

    def to_inchikey(self) -> str:
        """Standard InChIKey. (since 0.4)"""

    def to_json(self) -> str:
        """The .penz document (JSON). (since 0.4)"""

    def to_svg(self) -> str:
        """
        SVG, as the app exports it (the drawing embedded, so it reopens editable). (since 0.4)
        """

    def to_png(self, dpi: float = 300) -> bytes:
        """PNG bytes at this resolution. (since 0.4)"""

    @property
    def formula(self) -> str:
        """Molecular formula, Hill order, all molecules together (since 0.4)"""

    @property
    def mw(self) -> float:
        """Molecular weight (since 0.4)"""

    @property
    def exact_mass(self) -> float:
        """Monoisotopic mass (since 0.4)"""

    def save(self, path: str) -> None:
        """
        Write .penz (full fidelity), .mol, or ChemDraw .cdxml or .cdx; images go through export(). (since 0.4)
        """

    def export(self, path: str) -> None:
        """Write .svg, .png or .pdf, exactly as the app exports. (since 0.4)"""

    def __repr__(self) -> str: ...

def from_smiles(smiles: str) -> Document:
    """A structure from SMILES, laid out in 2D. (since 0.4)"""

def from_json(json: str) -> Document:
    """A Document from .penz JSON (as to_json writes). (since 0.4)"""

def read(path: str) -> Document:
    """
    Open a file: .penz, MOL, ChemDraw .cdxml/.cdx, .rxn, a Penzene SVG/PNG, or every record of an SDF, .smi or .inchi file laid out as a grid; records that can't be read are left out with a UserWarning naming them. (since 0.4)
    """

def drawing_styles() -> list[str]:
    """The drawing style names Document.style takes. (since 0.4)"""
