# Figures and exports

## Drawing styles

**Structure → Drawing Style** sets a document's style. Exports
come out at the style's real size. Values come from each publisher's ChemDraw stationery
or author instructions.

| Style | Bond length | Line | Bold | Double-bond gap | Font |
|---|---|---|---|---|---|
| ACS 1996 (default) | 14.4 pt | 0.6 pt | 2.0 pt | 18% | Arial 10 pt |
| JDP | 14.17 pt (0.5 cm) | 0.879 pt | 1.814 pt | 18% | IBM Plex Sans Light 10 pt |
| RSC | 12.2 pt | 0.449 pt | 1.602 pt | 20% | Helvetica 7 pt |
| Nature | 10.8 pt | 0.595 pt | 1.559 pt | 18% | Helvetica 6 pt |
| Angewandte | 17 pt | 0.75 pt | 2.6 pt | 18% | Arial 12 pt |
| SYNLETT / SYNTHESIS | 10.12 pt (0.357 cm) | 0.567 pt | 1.389 pt | 18% | Arial 6 pt |

Preferences sets the style for new documents.

## Export, copy and print

- **File → Export** writes SVG, PNG or PDF of the selection, or of everything. **Copy** puts the same
  picture on the clipboard (PNG and SVG) together with MOL, SMILES and the Penzene drawing, so
  pasting into Word, PowerPoint, Keynote or ChemDraw each gets what it understands.
- Exported SVG and PNG files carry the drawing inside them. Open or paste one back into Penzene and
  it's editable again, not a picture.
- **Preferences** sets the PNG resolution, a clear or white background, a scale (e.g. 85% to fit a
  journal column) and a margin.
- **File → Print** prints at the same size as an export, centred on the page.
- The canvas can follow the system's light or dark mode or use a theme (View → Theme). Exports are
  always black on clear or white. **Export Theme** saves the current one as a small JSON file (paper,
  ink, accent, hotspot, error and window colours) to edit or share; **Import Theme** adds one to
  the menu, and it is remembered.

## Page mode

**View → Page** shows a page on the canvas: A4, US Letter, or an ACS or RSC single or double
column at the journal's maximum figure height. Lay a scheme out at its final size; Export, Copy and
Print then take the whole page. A selection still exports just itself.

## Reaction schemes

Draw molecules either side of a reaction arrow, with reagents over or under it, and use
**Structure → Arrange Scheme** to line them up. **Edit → Copy as Reaction SMILES** and saving as
`.rxn` treat everything before the arrow as reactants, things over or under it as agents, and
things after it as products.
