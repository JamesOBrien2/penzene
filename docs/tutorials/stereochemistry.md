# Show calculated stereo labels

Penzene can calculate **R/S** labels for stereocentres and **E/Z** labels for suitable double
bonds. First the drawing must say which arrangement you mean: use a wedge or hashed bond, a
recognised Fischer projection, or an imported structure with stereochemistry. A flat ordinary
structure does not tell Penzene which enantiomer to choose.

## Try a structure with known stereochemistry

1. Download {download}`D-glyceraldehyde <files/fischer-glyceraldehyde.penz>` and open it with
   **File → Open**. Its OH is on the right of the Fischer cross.
2. Choose **Format → Stereo Labels**. A check mark appears in the menu, and **(R)** appears beside
   the stereocentre. This display choice is saved with the drawing and appears in exports.
3. Open {download}`the mirror image <files/fischer-glyceraldehyde-mirror.penz>` and turn on
   **Format → Stereo Labels** there too. Its OH is on the left and its label is **(S)**. Both
   drawings have the same formula, **C₃H₆O₃**.

```{image} ../_static/tutorials/fischer-glyceraldehyde.svg
:alt: D-glyceraldehyde Fischer projection with the OH group on the right
:width: 200px
```

## Label a structure you drew

For an ordinary skeletal drawing, choose **Wedge** or **Hash** from the **Bonds** tools and put it
on a bond leaving the stereocentre. A solid wedge points towards you; a hashed bond points away.
Then turn on **Format → Stereo Labels**. Penzene calculates the label from the drawing rather than
asking you to type R or S. Changing the wedge direction can change the assignment.

If no label appears, check that the centre has four different substituents and that its spatial
arrangement is specified. **Structure → Check Structure** reports stereocentres that need a wedge
and wedges placed on atoms that are not stereocentres. For double bonds, E/Z requires two
different substituents at each end and an unambiguous drawing.

The menu action shows or hides the calculated labels; it does not select an enantiomer for you.
For more on drawing wedges and hashes, see [Drawing](../drawing.md#bonds). For a longer Fischer
exercise, see [Newman and Fischer projections](projections.md).
