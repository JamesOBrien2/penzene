# A publication figure

Turn a structure into a figure that's ready for a journal: the journal's drawing style, colour where
it helps, and an export at the size it will be printed.

```{image} ../_static/tutorials/caffeine-figure.svg
:alt: Caffeine in the RSC style with nitrogens in blue and oxygens in red
```

Example file: {download}`caffeine-figure.penz <files/caffeine-figure.penz>` (the finished figure).

## 1. The journal's style

Choose **Format → Drawing Style** and pick the journal's: **ACS 1996** (ACS journals and most
others), **RSC**, or **JDP**. Bond length, line width, double-bond spacing and font all change to
match, and exports come out at that real size, so the figure needs no scaling to meet the guidelines.
**Edit → Preferences** sets the style for new documents.

## 2. Colour

Choose the colour tool from **Select**. Click its swatch to pick a colour: the CPK palette has the
usual element colours (blue for nitrogen, red for oxygen). Then click atoms, bonds, arrows or text to
paint them, or select several and paint them all at once. An atom's label takes its colour. The
fill tool in **Rings** shades a ring.

Use colour to point at something, such as the atoms a reaction changes; black stays clearest for
the rest, and every colour should still read in greyscale.

## 3. Size it on the page

**Page → Page Size** shows a page on the canvas: an ACS or RSC single or double column, A4 or
Letter. Lay the figure out at its final width. Export then takes the whole page, so the figure has
the column's width with the journal's margins.

## 4. Export

**File → Export…** as PDF or SVG for vector output, the usual choice for journals and slides. For
PNG, set the resolution in **Edit → Preferences** (600 dpi is a common requirement for line art),
along with a clear or white background and a scale. Before you send it:

- [ ] the style is the journal's, with nothing scaled after export
- [ ] it fits a column, at the width it will print
- [ ] labels are legible at that size, and colours still read in greyscale
- [ ] **Structure → Check Structure** finds no valence or stereo problems

The command line applies a style while exporting, for a set of figures that should match:

```sh
penzene --render caffeine-figure.penz caffeine.pdf --drawing-style RSC
```

More: [Figures and exports](../figures.md).
