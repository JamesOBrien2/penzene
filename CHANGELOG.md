# Changelog

What's new in each release. The app shows its own version's section once after an update:
a bullet with a **bold title** is a highlight (with an optional `<!-- icon: name -->` from
`resources/whatsnew/`); the rest are listed as text.

## Unreleased

- Rotate in 3D is a tool in the Select flyout: drag a molecule out of the page and it keeps its stereochemistry. The two Turn Over commands are removed (Shift+Alt+drag and the arrow keys still turn it).
- View → Show Implicit Hydrogens is removed: labels always show their implicit hydrogens (OH, NH₂).

## 1.5.0 (2026-09-28)

- The nmrshiftdb2 Database License for the built-in NMR shift table is now in the bundled third-party notices, so About → Third-party licenses shows it.
- **Predicted NMR shifts**: View → Predicted NMR Shifts writes each carbon's ¹³C shift, and the ¹H shift of atoms with hydrogens, beside the structure, looked up by HOSE code in data from nmrshiftdb2 (about 2.4 ppm ¹³C and 0.27 ppm ¹H mean error). <!-- icon: atom -->
- Cancelling the structure warning during Cut keeps the drawing, as it does for Copy.
- Copy pastes into Word and PowerPoint as a picture with a plain ⌘V or Ctrl+V (it no longer carries the SMILES as text; Copy as SMILES does), and the copied PNG keeps the drawing inside it on every platform, so pasting it back into Penzene gives the editable structure (not from Office, which re-renders pictures it copies).
- A reaction scheme that wraps onto a new row keeps the link between rows: the last product of one row is the next row's reactant when exported.
- Structure → Invert Stereochemistry turns every wedge into a hash and every hash into a wedge, giving the enantiomer without redrawing it (the selection's molecules, or the whole drawing).
- Copy as SMILES and to_smiles() keep stereo groups: a drawing with &1 or or1 centres is written as CXSMILES (`C[C@H](N)C(=O)O |&1:1|`), one without them as plain SMILES as before.
- **Stereo groups**: Tag a stereocentre abs, &1 or or1 from its right-click menu, as in ChemDraw; drawn beside the centre and kept in MOL V3000 and ChemDraw files, and read from CXSMILES. <!-- icon: atom -->
- Open reads MDL RD files (.rdf): the steps open as one scheme, left to right, with a step's product carrying on as the next step's reactant.
- On Windows, Copy also offers an Enhanced Metafile, so Word and PowerPoint paste a vector picture that stays sharp when scaled.
- Preferences can check structures before export and copy: unassigned stereocentres and valence errors are listed first, with Export Anyway or Cancel.
- The Python API reference gives the version that added each function, class, method and property, and the stability policy spells out that argument names, argument order and return types are stable while repr() strings and exact image bytes are not.
- The question icon in prompts such as "Save changes to this document?" shows on dark themes like Catppuccin Mocha, where it was black on macOS.
- The chair keys (9 and 0) fuse the chair onto the bond, pointing away from the ring, so it no longer overlaps the ring, even on a second bond of the same ring.
- **Isotope patterns**: View → Mass Spec Panel draws the isotope pattern of the selection, or the whole drawing, as a stick spectrum for [M], [M+H]⁺, [M+Na]⁺ or [M−H]⁻, with the m/z of the main peaks. <!-- icon: chart-bar -->

## 1.4.0 (2026-09-28)

- ChemDraw text opens with a ^ only where it keeps a charge raised: Cu²⁺ comes in as Cu2+, not Cu^2+.
- Typing a charged abbreviation such as N3- keeps the group and adds the charge (the azide anion), and the charge is drawn with the label.
- A ChemDraw reagent label without an id no longer deletes the molecules on opening.
- A lone abbreviation in a ChemDraw file (Boc on its own, say) opens as a labelled atom, not as text.

## 1.3.0 (2026-09-27)

- **Orbitals**: s, p, lobe and hybrid orbitals in Shapes, outlined, shaded or with a gradient, in any colour. Right-click to bring one to the front or send it behind the drawing; ChemDraw files keep them. <!-- icon: shapes -->
- **Isotopes**: Type 13C, 18OH or D on an atom (d makes deuterium); drawn ¹³C and kept through SMILES, MOL and ChemDraw files, with the formula (C[13C]H6O, CDCl3) and masses to match. <!-- icon: atom -->
- Charges in text: NH4+, Cu2+ and [Fe(CN)6]3- are set with the charge as a superscript; ^ marks one outright (SO4^2-). ChemDraw files keep counts and charges as sub- and superscripts.
- More abbreviations: COOH, Bu, iBu, OEt, NHBoc, NMe2, SO3H, TIPS, MOM and THP.
- A "no reaction" arrow, crossed through the middle, among the arrow tools; it round-trips with ChemDraw.
- ChemDraw export keeps abbreviations such as Boc and OMe as ChemDraw nicknames, which ChemDraw shows as labels and still knows the chemistry of, instead of drawing them out. ChemDraw nicknames also open as labels.
- Multi-step schemes export every step: Copy as Reaction SMILES gives a line per arrow, and Save As an RD file (.rdf) writes one Rxnfile per arrow. Before, only the first arrow was exported.

## 1.2.0 (2026-09-27)

- **Tool rail**: All the tools in one rail on the left, grouped as Select, Bonds, Rings, Atoms, Arrows and Shapes; a group's tools open beside it, and can be pinned open. The Draw, Chemistry and Figure switch is gone. <!-- icon: layout-sidebar -->
- **White and teal**: The page is white and the window a cool grey, in place of the pale yellow paper; the logo is the teal ring on a white tile, and dark mode is a deeper charcoal. <!-- icon: palette -->
- **Pages**: A drawing can hold several pages, as tabs along the bottom of the window, each with its own undo history; Edit → Move to Page moves a selection between them. View → Grid and View → Rulers measure the drawing at its final size. <!-- icon: layout-bottombar -->
- **SDF export**: Save As an SD file writes each molecule on the page as its own record, ready for a compound library. <!-- icon: database -->
- **Open from Finder and the file manager**: On macOS and Linux, .penz files open in Penzene with a double-click, and MOL, SDF, SMILES, Rxnfile and ChemDraw files offer it under Open With. <!-- icon: file-import -->
- A SMILES file's "SMILES Name" header row is skipped, not read as an unreadable structure.
- Erasing an atom also removes neighbours it leaves with no bonds, as erasing a bond already did, instead of leaving lone methanes behind.
- ChemDraw files keep their label size relative to their bonds, so crowded drawings with small labels don't overlap.
- An arrow with a huge bend, as a damaged file can hold, no longer exhausts memory when drawn.
- Atom numbers no longer sit on the same atom's lone pairs, radicals or δ label.
- Python: `Document.save` replaces a file whole, so a failed save leaves the old one intact, and it refuses image names (use `export()`) rather than writing MOL text into them.
- A ChemDraw text of size 0 no longer saves a drawing Penzene can't reopen, and drawings saved that way open again.
- An SDF file whose first structure has no title opens with that structure.
- Radicals survive export and copying to ChemDraw.
- Saving a drawing opened from an SDF or SMILES library asks where to save it instead of overwriting the library, and Save refuses image names rather than writing MOL text into them.
- Exports and copies no longer cut off atom numbers, lone pairs, δ labels or stereo labels at the edge, and a stereo label sits clear of its atom's number.
- Each running Penzene keeps its own crash-recovery copy, so opening a second one (as Windows does for each double-clicked file) no longer deletes the first one's, even when a later one reuses its process ID.
- View → Fit to Window (Ctrl+0) zooms to the selection when there is one.
- `penzene --render` with one output file draws every structure in an SDF or SMILES file, laid out as Open does, not just the last; a structure it can't read is reported, and the command fails.
- MOL export writes aromatic rings as the drawn single and double bonds; bond type 4 is only for queries, and some programs rejected it.
- Large drawings and pages can be scrolled to in full; the canvas no longer stops 5000 points from the origin.
- Clean keeps aromatic rings' double bonds where they were drawn, and keeps bond colours.
- MOL and ChemDraw files open with aromatic rings' double bonds where the file drew them.
- Drag a structure file onto the window to open it.
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
