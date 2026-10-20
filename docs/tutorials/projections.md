# Newman and Fischer projections

Insert a projection from the Templates panel, edit it, and check what it means when copied as
SMILES. Penzene treats Fischer projections as molecular structures and Newman projections as
schematic drawings; this difference matters when you export chemistry.

Example files:

- {download}`newman-butane.penz <files/newman-butane.penz>`: anti butane, viewed along its central bond.
- {download}`fischer-glyceraldehyde.penz <files/fischer-glyceraldehyde.penz>`: D-glyceraldehyde.
- {download}`fischer-glyceraldehyde-mirror.penz <files/fischer-glyceraldehyde-mirror.penz>`: its mirror image.

## 1. Insert a Newman projection

Start an empty document with **File → New**. Open **View → Templates**, expand **Projections**, and
click **Newman: butane, anti (staggered)**. The front carbon is at the centre; the circle represents
the back carbon. Three spokes run from the centre, and three from the circle's edge. The two
**CH₃** groups are 180° apart in this anti conformer.

```{image} ../_static/tutorials/newman-butane.svg
:alt: A staggered anti-butane Newman projection with the methyl groups opposite one another
:width: 300px
```

The inserted drawing is selected. Choose **Arrange → Group** (⌘G or Ctrl+G) to make its circle,
spokes and labels move together. With **Select**, drag it to a new position. To change a substituent,
choose the text tool from **Shapes** and click its label, or double-click the label with Select.
The text editor takes `CH3`, `Et` or any other text; formula digits become subscripts.

Changing a Newman label changes the illustration. The circle, spokes and text do not define atoms
or stereocentres, so this drawing has no molecular formula or SMILES. To exchange the molecule's
chemistry, draw a separate bonded structure with wedges and hashes, or import its stereochemical
SMILES. Keep the Newman drawing as the view of a chosen conformation.

## 2. Insert a Fischer projection

In another empty document, choose **View → Templates → Projections → Fischer: D-glyceraldehyde**.
The aldehyde is at the top, CH₂OH at the bottom, and OH on the right. In a Fischer projection,
horizontal bonds point towards you and vertical bonds point away.

```{image} ../_static/tutorials/fischer-glyceraldehyde.svg
:alt: D-glyceraldehyde as a Fischer projection with CHO above, CH2OH below, OH right and H left
:width: 220px
```

These are real atoms and bonds, with CHO and CH₂OH abbreviations. The status bar reads **C₃H₆O₃**.
Point at an atom and press **Enter** to edit its label, as for an ordinary molecule. Use Select to
move the entire projection. Keep each stereocentre's four bonds exactly horizontal and vertical;
tilting or moving an arm away from the cross can stop Fischer recognition.

## 3. Read the stereochemistry back

Choose **Edit → Copy As → SMILES**. Paste into a text editor: the glyceraldehyde string contains
`@` or `@@`, indicating the stereocentre. The exact choice of marker depends on the order of atoms
in the SMILES; `@` alone does not mean R or S. You can paste the string into a new Penzene document
to see the same stereochemistry in an ordinary structure drawing.

To make the mirror image, select the whole Fischer projection with **Edit → Select All**, then
choose **Arrange → Flip Horizontal**. OH moves left and H right, giving L-glyceraldehyde. Copy its
SMILES again: its stereochemistry changes, while its formula stays **C₃H₆O₃**. The mirrored example
file lets you compare the two without changing your original. **Format → Stereo Labels** shows
the assigned R/S label; D-glyceraldehyde is R and its mirror is S.

```{image} ../_static/tutorials/fischer-glyceraldehyde-mirror.svg
:alt: L-glyceraldehyde after a horizontal flip, with OH on the left and H on the right
:width: 220px
```

Penzene recognises an aligned four-bond Fischer cross when converting the drawing to chemistry.
A generic crossed sketch, or a Newman illustration made from shapes, does not gain stereochemistry
just because it looks like a projection.

## 4. Save and export

Save the working drawing as `.penz`. **File → Export…** writes SVG, PDF or PNG for teaching notes
and figures. For a Fischer projection, **Copy As → SMILES** or a MOL export also carries its molecular
stereochemistry. For the Newman illustration, use a figure format or `.penz` to retain the drawing.

More: [Templates and rings](../drawing.md#rings-and-templates),
[Stereochemistry](../chemistry.md), and [File formats](../files.md).
