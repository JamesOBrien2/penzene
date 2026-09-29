# File formats

| Format | Open | Save | Notes |
|---|---|---|---|
| Penzene `.penz` | ✓ | ✓ | everything: molecules, arrows, text, shapes, styles |
| MDL MOL / SDF | ✓ | ✓ | V2000 and V3000; a multi-record SDF opens as a grid, and Save As `.sdf` writes one record per molecule |
| SMILES `.smi`, InChI `.inchi` | ✓ | | one per line, as a grid; paste either as text too |
| MDL Rxnfile `.rxn` | ✓ | ✓ | and reaction SMILES by copy and paste, a line per step |
| MDL RD file `.rdf` | ✓ | ✓ | a multi-step scheme, one Rxnfile per arrow; opens as one scheme, left to right |
| ChemDraw XML `.cdxml` | ✓ | ✓ | molecules, abbreviations (as ChemDraw nicknames), arrows (curved too), text, lines, boxes, ovals, lone pairs |
| ChemDraw `.cdx` | ✓ | ✓ | as CDXML, on every platform |
| SVG, PNG, PDF | Penzene's own | export | the drawing rides along, so they reopen editable |

Open any of these with File → Open, or drag the file onto Penzene's window. On macOS and Linux,
installing Penzene makes it the app for `.penz` files and offers it under Open With for the
structure formats above.

## The .penz format

A `.penz` file is JSON: `{"format": "penzene", "version": 2, "pages": [...]}`, each page holding
its `name`, `atoms`, `bonds`, arrows, texts, fills, brackets and page size. The drawing settings
every page shares (`style`, `carbonLabels` and the like) sit at the top, beside `pages`; a page
whose setting differs keeps its own. It's written to be read back exactly. The
[JSON Schema](_static/penz.schema.json) describes every field. Here is the smallest drawing, a
single bond:

```json
{
    "format": "penzene",
    "version": 2,
    "pages": [
        {
            "name": "Page 1",
            "atoms": [{"x": 0, "y": 0, "z": 6}, {"x": 14.4, "y": 0, "z": 8}],
            "bonds": [{"a": 0, "b": 1, "order": 1}]
        }
    ]
}
```

Coordinates are points, x to the right and y down, with bonds 14.4 points long. Atoms are
referred to by their index in `atoms`. A field left out takes its default: carbon, no charge,
a single bond and so on. Two atoms share at most one bond.

Version 1, which Penzene 1.x writes, is one page's drawing with the format, version and settings
beside its atoms: `{"format": "penzene", "version": 1, "atoms": [...], "bonds": [...], ...}`. Its
other pages (Penzene 1.2 on) follow in its `pages` array, each a version 1 document of its own.
The clipboard and the drawing inside exported SVG, PNG and PDF files are still version 1.
Penzene 1.x can't open version 2, so Save As offers **Penzene 1 (single page)**: a
version 1 file, which those versions open (before 1.2, the first page only).

### Versions

Every Penzene release reads the files of every earlier one. Examples from each release live in
`tests/data/penz/`, and the tests open all of them.

- **Readers ignore fields they don't know.** A newer Penzene can add a field (say, a new arrow
  setting), and an older one still opens the file and draws what it understands.
- **New fields are optional**, and leaving one out means what files meant before it existed.
- **`version` changes only for a breaking change**: a field that changes meaning, or one that can't
  be left out. A reader refuses a version it doesn't know rather than guess.

## ChemDraw

CDXML import keeps what the file draws: molecules with their labels and abbreviations, arrows, text,
lines, boxes, ovals and lone-pair symbols. Pasting from ChemDraw works on macOS (and on Windows
where ChemDraw offers CDXML). The file's label size relative to its bond length comes with it, so
crowded drawings with small labels look as drawn; choosing a drawing style goes back to that style's
own. Penzene doesn't import ChemDraw's orbitals or TLC plates.
