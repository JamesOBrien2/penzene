# Figures and exports

## Drawing styles

**Format → Drawing Style** sets a document's style. Exports
come out at the style's real size.

| Style | Bond length | Line | Bold | Double-bond gap | Font |
|---|---|---|---|---|---|
| ACS 1996 (default) | 14.4 pt | 0.6 pt | 2.0 pt | 18% | Arial 10 pt |
| JDP | 14.17 pt (0.5 cm) | 0.879 pt | 1.814 pt | 18% | IBM Plex Sans Light 10 pt |
| RSC | 12.2 pt | 0.449 pt | 1.602 pt | 20% | Helvetica 7 pt |

Preferences sets the style for new documents.

## Export, copy and print

::::{container} features-grid

```{feature} file-export
:title: Export and copy
**File → Export** writes SVG, PNG or PDF of the selection, or of everything. **Copy** puts the same
picture on the clipboard (vector PDF, PNG and SVG) together with MOL, CDX and the Penzene
drawing, so Word, PowerPoint, Keynote and ChemDraw each get what they understand, and Office and
Keynote get sharp vectors rather than pixels. It carries no plain text, so a plain paste into Office
gives the picture; **Edit → Copy As → SMILES** (⌥⌘C) gives the text.
```

```{feature} file-type-pdf
:title: Figures that reopen
Exported SVG, PNG and PDF files carry the drawing inside them (a PDF also has the MOL file
attached). Open or paste one back into Penzene and it's editable again, not a picture.
```

```{feature} adjustments
:title: Resolution, background and scale
**Preferences** sets the PNG resolution, a clear or white background, a scale (e.g. 85% to fit a
journal column) and a margin.
```

```{feature} printer
:title: Print
**File → Print** prints at the same size as an export, centred on the page.
```

```{feature} moon
:title: Light and dark
The canvas can follow the system's light or dark mode or use a theme (**View → Theme**). Exports
are always black on clear or white.
```

::::

## Editing a drawing in Word or PowerPoint (Windows)

Penzene installed with the Windows installer registers a **Penzene Drawing Object**, as ChemDraw
does. Copy a drawing and paste it into Word or PowerPoint: Ctrl+V gives the object, and
**Paste Special** also lists it by that name. You can also add a new, empty one with **Insert → Object**.
The document shows the picture and keeps the drawing (every page, as a `.penz`) inside the
`.docx` or `.pptx`, so the picture still shows on a computer without Penzene.

Double-click the object to open it in a Penzene window titled "Drawing in …". Edit it, then
**File → Save** (Ctrl+S) or close the window and choose Save: the picture in the document updates.
Save the Word or PowerPoint file to keep the change. The picture is page 1, drawn with the export
settings in Preferences.

The portable `.zip` doesn't register Penzene, so from it Copy pastes a picture. macOS Office has no
embedded objects; there, a pasted PDF or PNG from Penzene opens back up in Penzene (see above).

## Page mode

**Page → Page Size** shows a page on the canvas: A4, US Letter, or an ACS or RSC single or double
column at the journal's maximum figure height. Lay a scheme out at its final size; Export, Copy and
Print then take the whole page. A selection still exports just itself.

## Reaction schemes

Draw molecules either side of a reaction arrow, with reagents over or under it, and use
**Arrange → Arrange Scheme** to line them up. **Edit → Copy As → Reaction SMILES** and saving as
`.rxn` treat everything before the arrow as reactants, things over or under it as agents, and
things after it as products.

**Arrange → Number Compounds** puts a bold number under each selected molecule (or every molecule)
that has none. The numbers follow their molecules and stay in scheme order, row by row and left to
right, as you edit: insert a step and the later numbers move up. Edit one to add a letter; **2a** and
**2b** keep sharing their number. They run on across the pages in tab order, so a scheme split over
pages keeps one sequence.

```{image} _static/scheme-light.png
:alt: A reaction scheme: salicylic acid to aspirin with Ac₂O and H₂SO₄ over the arrow
:class: shot only-light
```

```{image} _static/scheme-dark.png
:alt: A reaction scheme: salicylic acid to aspirin with Ac₂O and H₂SO₄ over the arrow, dark theme
:class: shot only-dark
```
