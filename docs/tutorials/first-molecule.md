# Your first molecule in five minutes

Draw ethanol, save an editable drawing, and export a figure. Once Penzene is installed, this takes
about five minutes. You need no chemistry file or online lookup.

Example file: {download}`first-molecule.penz <files/first-molecule.penz>` (the finished drawing).

## 1. Install and start

Download the installer for your computer from the
[latest release](https://github.com/JamesOBrien2/penzene/releases/latest) and follow
[Installing](../install.md). Open Penzene, then choose **File → New** for an empty drawing.
If a drawing is already open, save it when prompted.

## 2. Draw two bonds

Click **Bonds** in the rail on the left and choose the **Single bond** tool (the plain line).
Alternatively, press Esc, then `x`.

1. Click empty space near the centre of the canvas. A single bond appears: its two ends are carbons.
2. Click the right-hand end of that bond. A second bond grows from it, making a three-atom zigzag.
3. Point at the new end until it lights up, then type lowercase `o`. That carbon becomes oxygen,
   displayed as **OH**. Do not click it with the bond tool, which would add another bond.

The two unlabelled positions are carbons; Penzene supplies their hydrogens automatically. The status
bar should read **C₂H₆O**, with a molecular weight of about **46.07**. If it does not, use **Edit → Undo**
(⌘Z on macOS, Ctrl+Z elsewhere) and repeat the last step. Press Esc to clear the hotspot.

```{image} ../_static/tutorials/first-molecule.png
:alt: Penzene showing the finished two-bond ethanol drawing and C2H6O in the status bar
:width: 800px
```

## 3. Save the drawing

Choose **File → Save** (⌘S or Ctrl+S), choose the **Penzene** file type, and name the file
`ethanol.penz`. This keeps the atoms, bonds and drawing settings editable. Open the file with
**File → Open** whenever you want to continue.

## 4. Export a figure

Choose **File → Export…** (⌘E or Ctrl+E), select **SVG**, and save `ethanol.svg`. SVG keeps the
lines and letters sharp when resized. PDF is another vector option; choose PNG when you need a
bitmap, with its resolution and background set in **Edit → Preferences**.

Here is the exported figure, enlarged for display:

```{image} ../_static/tutorials/first-molecule.svg
:alt: Ethanol exported from Penzene as a sharp vector drawing
:width: 180px
```

The exported SVG also reopens in Penzene as editable atoms and bonds. Keep the `.penz` file as your
working document, and use the exported figure in your slides or manuscript.

Next, try [Caffeine from the keyboard](keyboard.md) or [A reaction scheme](reaction-scheme.md).
For all the tools and hotkeys, see [Drawing](../drawing.md) and [Keyboard shortcuts](../keys.md).
