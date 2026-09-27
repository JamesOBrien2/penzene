# Changelog

What's new in each release. The app shows its own version's section once after an update:
a bullet with a **bold title** is a highlight (with an optional `<!-- icon: name -->` from
`resources/whatsnew/`); the rest are listed as text.

## Unreleased

- **Tool rail**: All the tools in one rail on the left, grouped as Select, Bonds, Rings, Atoms, Arrows and Shapes; a group's tools open beside it, and can be pinned open. The Draw, Chemistry and Figure switch is gone. <!-- icon: layout-sidebar -->
- **White and teal**: The page is white and the window a cool grey, in place of the pale yellow paper; the logo is the teal ring on a white tile, and dark mode is a deeper charcoal. <!-- icon: palette -->
- **Pages**: A drawing can hold several pages, as tabs along the bottom of the window, each with its own undo history; Edit → Move to Page moves a selection between them. View → Grid and View → Rulers measure the drawing at its final size. <!-- icon: layout-bottombar -->
- **SDF export**: Save As an SD file writes each molecule on the page as its own record, ready for a compound library. <!-- icon: database -->
- ChemDraw files keep their label size relative to their bonds, so crowded drawings with small labels don't overlap.
- An arrow with a huge bend, as a damaged file can hold, no longer exhausts memory when drawn.
- Python: `Document.save` replaces a file whole, so a failed save leaves the old one intact, and it refuses image names (use `export()`) rather than writing MOL text into them.
- A ChemDraw text of size 0 no longer saves a drawing Penzene can't reopen, and drawings saved that way open again.
- **Open from Finder and the file manager**: On macOS and Linux, .penz files open in Penzene with a double-click, and MOL, SDF, SMILES, Rxnfile and ChemDraw files offer it under Open With. <!-- icon: file-import -->
- An SDF file whose first structure has no title opens with that structure.
- Radicals survive export and copying to ChemDraw.
- Saving a drawing opened from an SDF or SMILES library asks where to save it instead of overwriting the library, and Save refuses image names rather than writing MOL text into them.
- Each running Penzene keeps its own crash-recovery copy, so opening a second one (as Windows does for each double-clicked file) no longer deletes the first one's.
- `penzene --render` with one output file draws every structure in an SDF or SMILES file, laid out as Open does, not just the last.
- MOL export writes aromatic rings as the drawn single and double bonds; bond type 4 is only for queries, and some programs rejected it.
- Clean keeps aromatic rings' double bonds where they were drawn, and keeps bond colours.
- MOL and ChemDraw files open with aromatic rings' double bonds where the file drew them.
- Typing a charged label such as NH3+, O- or Fe3+ sets the element and its charge, not an unknown group.

## 1.1.1 (2026-09-27)

- A .penz file with an unknown element or two bonds between the same atoms is refused, instead of closing Penzene.
- After an erase or an undo, the selection no longer lands on other atoms (where Delete would remove them); opening a file starts with nothing selected.
- The atom tool replaces an abbreviation: clicking Ph with N gives N, where before the Ph stayed.
- A click on a selection that wobbles a pixel no longer nudges it (without an undo step).
- Commands with nothing to do (Remove Explicit Hydrogens with none, Clean on a clean drawing) no longer add an undo step or mark the file unsaved.
- Saving writes a complete copy before replacing the file, so a save that fails part way (a full disk) no longer leaves it cut short; the same goes for autosave.
- A charged formula shows its charge as a superscript (O₄S²⁻ for sulfate), in the status bar and the Properties panel.
- Delete on an attachment point removes it, where before it became a methyl.
- A dashed arrow's head is drawn solid, not with a broken outline.
- Dragging an atom onto another keeps its brackets and its ring's circle setting, where before they lost it.
- Image exports (PNG, SVG, PDF) and descriptor tables replace a file only once they are written in full, so a failed export no longer leaves it cut short.

## 1.1.0 (2026-09-26)

- **Descriptor tables**: File → Export Descriptors… and penzene --descriptors write a CSV with one row per molecule (SMILES, InChIKey, mass, cLogP, TPSA, H-bond donors and acceptors, and more), keeping unreadable structures as rows that say why. <!-- icon: table -->

## 1.0.0 (2026-09-26)

- The .penz format is documented, with a JSON Schema (penzene.readthedocs.io, File formats), and every release's files are tested to keep opening.
- The j hotkey (η⁵-cyclopentadienyl) on a ring that's already there makes that ring Cp⁻, instead of corrupting memory; Python scripts no longer crash now and then as they exit.
- The Drawing, Chemistry and Figures pages of the documentation show the app, in light and dark.
- Template thumbnails are drawn in the theme's ink, so they show on dark panels, and they're sharper.
- Bold bonds join without a notch where they meet other bonds, on screen and in exports.
- Copy as SMILES gives nothing, rather than SMILES other programs can't read, for impossible structures such as a hydrogen with four bonds.

## 0.9.0 (2026-09-26)

- **Lab notebook look**: warm paper and teal controls, with tools grouped into Draw, Chemistry and Figure. <!-- icon: palette -->
- **Welcome**: an empty page offers examples (aspirin, a reaction scheme, a mechanism) and links to the keys and documentation. <!-- icon: sparkles -->
- **Documentation**: a full guide online at penzene.readthedocs.io, in light and dark. <!-- icon: book-2 -->
- **Updates**: Help → Check for Updates, and an optional weekly check (off unless you turn it on). <!-- icon: refresh -->
- **Sharp copies**: Copy puts a vector PDF on the clipboard for Word, PowerPoint and Keynote, and every PDF carries the drawing, so it opens or pastes back editable. <!-- icon: file-type-pdf -->
- **ChemDraw files everywhere**: binary .cdx opens and saves on Windows too, and keeps arrows and text. <!-- icon: file-import -->
- The logo is drawn in the same teal on warm paper as the app.
- Large drawings stay quick: drawing and chemistry now scale linearly (a 2000-atom page edits in about 17 ms).
- Preferences → Language, and the groundwork for translations (see Contributing to add one).
- Accessibility: Tab reaches every tool, with a visible focus ring and names for screen readers; arrow keys move around the periodic table; the hotspot is announced; every theme meets WCAG contrast.
- The What's New window highlights main additions with icons and lists smaller changes below.
- ChemDraw files: plain lines are no longer imported twice; lone pairs come through.
- The colour tool offers the CPK atom colours; click its swatch to pick one, then click atoms and bonds.
- Python: type stubs for editors and type checkers, and a written stability policy (semantic versioning, one minor release of deprecation warnings before anything is removed).
- Python: save() writes ChemDraw .cdxml and .cdx too.
- Arrow heads sit evenly on tightly curved arrows, and half heads (fishhook, equilibrium) are clean at the base.
- A dashed ellipse joins the Figure tools, which now pair each solid shape with its dashed version.
- A round handle above a selection rotates it; hold Shift for 15° steps, or Ctrl to land it square to the page or at 45°.
- The downloads build RDKit without its ChemDraw library (which contains MPL-1.0 code), so everything they ship is GPL-compatible; Penzene reads and writes ChemDraw files with its own copy.

## 0.8.0 (2026-09-25)

- **Template library**: amino acids, sugars, nucleobases, scaffolds, and your own (View → Templates). <!-- icon: books -->
- **Projections**: Haworth, Fischer and Newman drawings; Haworth and Fischer give the right stereo. <!-- icon: hexagons -->
- **Rotate in 3D**: Shift+Alt+drag turns a structure out of the page and keeps its stereo. <!-- icon: rotate-3d -->
- **Interactions and transition states**: dotted H-bonds and dashed forming or breaking bonds. <!-- icon: line-dashed -->
- **Arrange Scheme**: lines a reaction scheme up, with reagents centred over their arrows. <!-- icon: layout-distribute-horizontal -->
- Lone pairs, radicals, δ+ and δ−, and brackets with a subscript.
- Lines, boxes, rounded boxes and ellipses.
- Stretch and squash with handles; Structure → Transform.
- R-groups (R1, R2…) and generic atoms (X, Ar); attachment points and η-bonded rings.
- Atom Properties; bring a bond to the front so crossings show a gap.
- Arrow keys nudge a selection; Space and Enter move between the hotspot and its molecule.

## 0.7.0 (2026-09-24)

- Exported SVG and PNG reopen as editable drawings.
- Reactions: copy as reaction SMILES, open and save .rxn.
- Multi-record SDF, .smi and .inchi files open as a grid; MOL V3000; paste InChI.
- Save as ChemDraw CDXML (and CDX); paste from ChemDraw.
- Export scale and margin; page mode (A4, Letter, journal columns); printing.

## 0.6.0 (2026-09-24)

- Stereo labels: CIP (R)/(S) and (E)/(Z).
- Check Structure explains valence, stereo and label problems.
- Properties panel: cLogP, TPSA, H-bond donors and acceptors, rotatable bonds, Lipinski, Veber.
- Carbon and hydrogen display options; add or remove explicit hydrogens.
- Aromatic circles; atom numbers and reaction atom-map numbers.
- X and R on an atom, and free-text atom labels.
- Name to structure and structure to name via PubChem (online, on request).
- The version in the window title; `pip install penzene` from PyPI.

## 0.5.0 (2026-09-24)

- Preferences; recent files, autosave and crash recovery.
- Right-click menus; flip, align and distribute; colour atoms, bonds and text.
- Drag an atom onto another to merge; drag to size rings.
- Per-style bond lengths; fixes to Clean and to CDXML import.

## 0.4.0 (2026-09-24)

- The Python package: `import penzene`, with wheels for macOS, Linux and Windows.
- `penzene --render` for batch figures from the command line.

## 0.3.0 (2026-09-24)

- ACS, RSC and JDP drawing styles; ring fill; dark mode and Catppuccin themes.
- A new tool palette with a periodic-table picker.

## 0.2.0 (2026-09-24)

- Reaction arrows and text; abbreviations (Me, Ph, Boc…).
- Formula and molecular weight; InChI and InChIKey; CDXML import.
- Windows installer and an Intel Mac build.

## 0.1.0 (2026-09-23)

- The first release: draw with bonds, chains, rings and keyboard shortcuts; undo and redo.
- MOL, SDF, SMILES and .penz; Clean; implicit hydrogens.
- SVG, PNG and PDF export; copy as image, MOL and SMILES.
