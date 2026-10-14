# Journal-ready export

Prepare a chemical figure at its intended physical size, choose vector or raster output, and check
the file you will submit. Start with the journal's current author instructions: its required
dimensions, fonts and resolution can differ from the presets.

Example files:

- {download}`journal-caffeine.penz <files/journal-caffeine.penz>`: black caffeine in RSC style, without a page.
- {download}`journal-caffeine.pdf <files/journal-caffeine.pdf>` and
  {download}`journal-caffeine-600dpi.png <files/journal-caffeine-600dpi.png>`: its cropped exports.
- {download}`journal-caffeine-column.penz <files/journal-caffeine-column.penz>` and
  {download}`journal-caffeine-column.pdf <files/journal-caffeine-column.pdf>`: the same drawing on an RSC single-column page.

```{image} ../_static/tutorials/journal-caffeine.svg
:alt: A black caffeine structure in the RSC drawing style, tightly cropped for export
:width: 320px
```

The preview above is enlarged for reading; its displayed size in this web page is not its print size.

## 1. Choose a style before laying out the figure

Open `journal-caffeine.penz`. Choose **Format → Drawing Style → RSC** to reproduce the example.
Styles set physical bond length, line widths, double-bond spacing and the label font together.

| Built-in style | Bond length | Label font setting |
|---|---|---|
| ACS 1996 | 14.4 pt (5.08 mm) | Arial, 10 pt |
| RSC | 12.2 pt (about 4.30 mm) | Helvetica, 7 pt |
| JDP | 14.17 pt (about 5.00 mm) | IBM Plex Sans Light, 10 pt |

These are the preset settings, rather than a promise of compliance with every journal's current
rules. Check the applicable instructions. **Edit → Preferences** chooses the style for new
documents; **Format → Drawing Style** changes this document.

Fonts must be available on the computer generating the figure. A substitute font can change label
widths and spacing. Penzene exports letters as outlines in SVG and PDF, so their appearance is kept
on the receiving computer, but that cannot correct a substitution that happened before export.
Inspect the generated file, including Greek letters, subscripts and charges.

## 2. Choose the exported area

With **Page → Page Size → None**, export crops to the drawing. This is useful for a single structure
or a compact scheme, as in `journal-caffeine.pdf`. The cropped width depends on the figure's content.
**Edit → Preferences → Margin around exports** adds a margin in points around a cropped figure.

For a figure laid out within a column, choose **Page → Page Size → RSC single column**, or an ACS/RSC
double-column preset. Export takes the entire page when nothing is selected. The RSC single-column
preset is approximately **83 mm wide and 233 mm high**; the large blank area in the column example
comes from exporting that full-height page. Page mode uses the page itself as the frame, rather than
adding the crop margin.

Press Esc before exporting the full figure. Exporting a selection crops to its contents, even if
the document has a page. Keep the drawing inside the page boundary: material beyond it is clipped.
If you need a compact figure instead of the full-height page, choose **None** and export the crop.

## 3. Keep scale and resolution separate

In **Edit → Preferences**, set **Export and copy scale** to **100%** to preserve the style and page
dimensions. Reducing the scale also reduces bond lengths, line widths, fonts and the output page.
Zooming the canvas only changes your view; it does not change the exported size.

For PNG, set **PNG resolution** to the journal's requirement. At a fixed print size, higher DPI gives
more pixels, without making the printed labels larger. The relation is:

```text
pixels = width in millimetres / 25.4 × DPI
```

An 83 mm-wide figure needs about **1961 pixels at 600 dpi**, or **980 pixels at 300 dpi**, with
rounding to whole pixels. The column example was checked at 600 dpi and is 1961 pixels wide.
Set **Export and copy background** to **White** when the destination needs an opaque figure;
**Clear** gives transparency. A screenshot of the canvas does not have the export's calibrated size.

## 4. Export and inspect the actual file

Choose **File → Export…** and save PDF or SVG when the journal accepts vector figures. Their lines
stay sharp when displayed or printed at another resolution, but resizing them still changes their
physical line widths and label sizes. Set the intended insertion size in your manuscript and keep
the aspect ratio.

For PNG, choose the resolution first in Preferences, then export. Check its pixel dimensions and
print dimensions in an image viewer. Exported dimensions can differ slightly from the nominal
size because of rounding. For PDF, check the page size in a PDF viewer and print at **Actual size**
or **100%**, rather than fitting it to another page.

For comparison, the column example's PDF records **235 × 661 pt** (about **82.9 × 233.2 mm**),
while its 600 dpi PNG is **1961 × 5504 pixels**. Check the dimensions of the actual output rather
than relying only on the preset's name.

The example's PDF and PNG were exported by Penzene's rendering engine. The PDF is vector output;
the PNG records approximately 600 dpi. Open either back in Penzene to recover the editable drawing.
Keep the `.penz` source alongside the submitted figure.

## 5. Before submission

- [ ] The style, font and physical dimensions agree with the journal's current instructions.
- [ ] **Structure → Check Structure** reports no unresolved valence, overlap or stereo problems.
- [ ] Wedges, hashes, charges, abbreviations and reaction conditions say what you intend.
- [ ] The exported area includes the whole figure, with an appropriate crop or page frame.
- [ ] Labels remain legible at the final insertion size; colour distinctions also work in greyscale.
- [ ] The PNG has the required pixels and DPI, or the vector file has the intended page size.
- [ ] The destination has not stretched the figure or resized its labels below the required size.
- [ ] You have reopened the final export and kept its editable `.penz` source.

For a batch, the command line reproduces the style without a window:

```sh
penzene --render journal-caffeine.penz journal-caffeine.pdf --drawing-style RSC
```

More: [A publication figure](publication-figure.md) and [Figures and exports](../figures.md).
