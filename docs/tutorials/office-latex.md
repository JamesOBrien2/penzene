# Figures in Word, PowerPoint and LaTeX

Put the same caffeine figure into a manuscript, a slide and a LaTeX document. Keep its
`.penz` source beside the manuscript so you can revise the chemistry later.

Download these examples into one folder:

- {download}`caffeine-figure.penz <files/caffeine-figure.penz>`: the editable source from
  [A publication figure](publication-figure.md), with the RSC style and coloured N and O atoms.
- {download}`office-caffeine.svg <files/office-caffeine.svg>`: a vector picture for Word and PowerPoint.
- {download}`office-caffeine.pdf <files/office-caffeine.pdf>`: the same figure for LaTeX.
- {download}`office-figure.tex <files/office-figure.tex>`: a complete LaTeX document using that PDF.

```{image} files/office-caffeine.svg
:alt: Caffeine with blue nitrogen labels and red oxygen labels, exported in the RSC drawing style
```

## 1. Prepare the figure

Open `caffeine-figure.penz` in Penzene. Click empty canvas to clear any selection: **Copy** and
**Export** use the selection when there is one, otherwise the current page. Choose
**Structure → Check Structure** and keep the source with **File → Save As…**.

Choose **File → Export…** and save `office-caffeine.svg` and `office-caffeine.pdf`. Both are
vectors: enlarging the picture does not reveal pixels. The exports use the drawing style's
physical bond length and label size, independently of the canvas zoom. The downloads were
exported from the actual Penzene renderer, rather than captured from the screen.

For this small figure, leave **Page → Page Size** off so the export fits the drawing. Page mode
exports the whole page, including its empty space. See [A publication figure](publication-figure.md)
for journal sizes, colours and margins.

## 2. Insert a vector picture in Word or PowerPoint

In a recent desktop version of Word or PowerPoint, choose **Insert → Pictures** and the option
for a file on your computer, then select `office-caffeine.svg`. In Word, **In Line with Text**
keeps a small structure with its paragraph; add a caption through **References → Insert Caption**.
In PowerPoint, place the picture beside your explanation and add the caption in a text box.

Keep the aspect ratio locked when resizing. Check the labels at the size the reader will see.
Resizing also changes bond lengths and font sizes, so keep the export's original size when a
journal requires the RSC dimensions. Change the drawing in Penzene and export again if it needs
more space; stretching the picture changes its proportions.

Office versions that cannot insert SVG can use a PNG exported from Penzene. Set its resolution
in **Edit → Preferences** before exporting; 600 dpi is a useful starting point for chemical line
art. PNG is a bitmap, so check it at its final size. Use SVG for supported Office versions and PDF
for the LaTeX workflow below, rather than assuming Office will insert a PDF as a picture.

## 3. Copy and paste directly

For a quick transfer, clear the selection or select just the figure, then choose **Edit → Copy**
(Ctrl+C on Windows/Linux, Cmd+C on macOS). In Word or PowerPoint, paste with Ctrl+V or Cmd+V.
Use ordinary Copy, rather than **Copy As → SMILES**, which copies chemical text.

Penzene offers several clipboard representations; the receiving app chooses one. On macOS it
offers PDF and SVG vectors as well as PNG. On Windows it also offers an EMF vector picture and,
when installed with the installer, the editable drawing object described below. **Paste Special**
lets you choose a picture instead of an object where Office offers both. File insertion from step
2 gives you explicit control over the format if a paste arrives as a bitmap.

Save the `.docx` or `.pptx`, close it and reopen it to check the figure. Keep the `.penz` source:
copying a picture back out of Office may lose the embedded chemistry. In particular, PowerPoint
on macOS re-renders the picture and drops that drawing data. Reopening the original Penzene SVG
or PDF is editable; copying an Office picture back is not a reliable replacement for the source.
To revise a static picture, edit the `.penz`, export again and replace the picture in Office.

## 4. Keep an editable drawing object on Windows

Use desktop Word or PowerPoint on Windows with Penzene installed by the **Windows installer**.
The portable ZIP does not register the **Penzene Drawing Object**. Office for Mac and Office on
the web do not provide this Penzene OLE editing workflow.

1. Open the sample in Penzene, clear the selection and **Copy** (Ctrl+C).
2. In Word or PowerPoint, **Paste** (Ctrl+V), or use **Paste Special → Penzene Drawing Object**.
   Choose Paste, not Paste Link: this keeps a copy of the drawing inside the Office file.
3. Double-click the object. Penzene opens a separate window titled **Drawing in …**. For a visible
   trial edit, change one atom's colour using the Colour tool in the Select group.
4. Choose **File → Save** (Ctrl+S) in that Penzene window. The picture in Office updates. You can
   also close the editing window and choose Save when prompted.
5. Save the Word or PowerPoint document too. Close and reopen it, then double-click the object
   to check that the drawing still contains your edit.

**Insert → Object → Penzene Drawing Object** starts an empty embedded drawing instead. An
embedded drawing can contain multiple pages; Office displays page 1. Copy initially takes the
current page or selection, not all the tabs in your source document. A computer without Penzene
can show the cached picture, but needs Penzene installed to edit its chemistry.

If Paste Special has no Penzene object, check that you used the installer and copied from
Penzene rather than from an exported image. A static picture cannot become an editable OLE
object just by double-clicking it. Use the vector-picture workflow if your Office setup does not
offer the object. This uses Penzene's existing Windows integration.

## 5. Include the PDF in LaTeX

Put `office-figure.tex` and `office-caffeine.pdf` in the same folder, or upload both to an Overleaf
project. Set `office-figure.tex` as the main document and use pdfLaTeX. Locally, run this from
that folder (the command works in fish too):

```sh
pdflatex -halt-on-error office-figure.tex
```

Run the command twice to resolve the figure reference. The example loads `graphicx` and
includes the PDF at its natural size:

```latex
\begin{figure}[htbp]
  \centering
  \includegraphics{office-caffeine.pdf}
  \caption{Caffeine drawn in Penzene using the RSC drawing style.}
  \label{fig:caffeine}
\end{figure}
```

PDF works directly with pdfLaTeX, LuaLaTeX and XeLaTeX; this example needs no SVG converter or
shell escape. If you set a width, for example `\includegraphics[width=35mm]{office-caffeine.pdf}`,
the aspect ratio stays intact, but the bond lengths and labels scale too. Open the compiled PDF
and check its final size. The caption belongs in LaTeX, so it can be numbered and referenced with
`Figure~\ref{fig:caffeine}`.

To update the figure, edit the saved `.penz`, export over `office-caffeine.pdf` and compile again.
Keep the source, exported figure and `.tex` together; the final manuscript PDF should not be
used as the editable chemistry archive.

## Try it

Change the oxygen labels to black in the source. Replace the SVG picture in Word or PowerPoint,
or save the change through a Windows embedded object's Penzene window. Export the PDF again
and recompile the LaTeX example. Each document should show the same change, with no stretched
bonds or clipped labels.

More: [Figures and exports](../figures.md) and [Files and interoperability](../files.md).
