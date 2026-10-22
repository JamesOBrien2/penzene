# Show calculated stereo labels

Two drawings can have the same atoms and bonds but differ in how those atoms are arranged in
space. Penzene can calculate **R/S** labels for stereocentres and **E/Z** labels for suitable double
bonds. **Format → Stereo Labels** shows or hides the results. It does not choose a configuration
for you: the drawing must first specify one.

## Read wedges, hashes and wavy bonds

| Bond | What it says |
|---|---|
| Plain line | The bond is drawn in the plane of the page. On its own, it does not specify which enantiomer you mean. |
| Solid wedge | The bond points **towards** you, out of the page. |
| Hashed wedge | The bond points **away** from you, behind the page. |
| Wavy bond | The bond's spatial direction is **unspecified**. Use it when you do not know which arrangement is present. |

For a tetrahedral stereocentre, the **narrow end** of a wedge or hash starts at the stereocentre.
The broad end reaches the group whose direction it shows. A hashed wedge is different from the
ordinary dashed bond tool, which only changes the line's appearance. Penzene's **Bonds** tools
include Wedge, Hash and Wavy; when your pointer is over an existing bond, `w`, `h` and `y` apply
those styles too.

## Draw a stereocentre yourself

Try 2-butanol. It has one stereocentre: its second carbon is attached to **OH**, **H**, **CH₃** and
**CH₂CH₃**. Those four groups are different.

1. Choose **File → New**, then **Bonds → Single bond**. Click empty space to make one bond. Click
   its right-hand end twice more to make a four-carbon zigzag.
2. Choose **Bonds → Wedge**. Start at the **second carbon from the left** and drag into empty space
   below the chain to add a branch. The narrow end should touch that carbon.
3. Point at the new branch end and type lowercase `o`. The atom becomes **OH**.
4. Choose **Format → Stereo Labels**. An **(R)** or **(S)** label appears beside the second carbon.
   The answer depends on the arrangement you drew; Penzene calculates it rather than asking you
   to type a letter.
5. Choose **Bonds → Hash** and click the carbon–OH bond. The OH now points behind the page instead
   of towards you, so the R/S label should switch. You can switch it back with Wedge.

Save the drawing as `.penz` if you want to keep it. The checked Stereo Labels setting is saved with
the document and the calculated labels appear in figure exports.

## Try a known example without a download

For a result you can compare with a named molecule, choose **File → New**, open **View →
Templates**, expand **Projections**, and click **Fischer: D-glyceraldehyde**. Its OH is on the right
of the cross. Choose **Format → Stereo Labels**: the stereocentre is **(R)**.

```{image} ../_static/tutorials/fischer-glyceraldehyde.svg
:alt: D-glyceraldehyde Fischer projection before stereo labels are switched on, with OH on the right
:width: 200px
```

In a Fischer projection, horizontal bonds point towards you and vertical bonds point away. Select
the whole drawing with **Edit → Select All**, then choose **Arrange → Flip Horizontal**. OH moves to
the left and the label becomes **(S)**. Both forms have formula **C₃H₆O₃**. Keep the four arms
horizontal and vertical: tilting the cross can stop Penzene recognising it as a Fischer projection.

If you prefer to start from files, the {download}`D-glyceraldehyde drawing
<files/fischer-glyceraldehyde.penz>` and its {download}`mirror image
<files/fischer-glyceraldehyde-mirror.penz>` are optional examples.

## What R/S and E/Z mean

**R/S** names the arrangement at a stereocentre using Cahn–Ingold–Prelog priorities. The atom
directly attached with the higher atomic number has higher priority; if two are tied, compare the
next atoms outward until they differ. With the lowest-priority group pointing away, priorities
1 → 2 → 3 run clockwise for **R** and anticlockwise for **S**. A wedge alone does not mean R, and a
hash alone does not mean S. R/S also does not tell you whether a substance rotates light to the
right or left.

**E/Z** describes a double bond when each end has two different groups. Choose the higher-priority
group at each end using the same priority rules. If those groups are on the **same** side, the bond
is **Z**; if they are on **opposite** sides, it is **E**. Penzene reads E/Z from the drawn double-bond
geometry, so you do not need to add a wedge to an alkene.

## If no label appears

Check that **Format → Stereo Labels** has a check mark, and that the centre really has four
different groups. A flat tetrahedral centre with no wedge or hash has no R/S configuration to
calculate. A double bond needs two different groups at **each** end and an unambiguous drawing.
Use **Structure → Check Structure** to find possible stereocentres without a wedge and wedges that
start at non-stereocentres.

For more drawing controls, see [Bonds](../drawing.md#bonds). For a longer projection exercise, see
[Newman and Fischer projections](projections.md).
