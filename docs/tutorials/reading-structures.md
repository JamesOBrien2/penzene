# Reading a structure drawing

Chemical drawings leave some atoms and hydrogens out to keep the picture clear. This short guide
explains the shorthand you need for the [first molecule tutorial](first-molecule.md). You do not
need to know SMILES or use a keyboard shortcut.

```{image} ../_static/tutorials/first-molecule.svg
:alt: Ethanol: two connected bonds ending in an OH label
:width: 220px
```

This is ethanol. Read it from left to right:

| What you see | What it means here |
|---|---|
| An unlabelled line end | A carbon atom |
| An unlabelled corner where the lines meet | Another carbon atom |
| `OH` at the right-hand end | An oxygen atom with one hydrogen |
| A line between two positions | A single bond between those atoms |

The left carbon has three hydrogens and the middle carbon has two. Penzene supplies those
hydrogens automatically, so the whole drawing is **CH₃–CH₂–OH**, or **C₂H₆O**. Two parallel lines
between atoms would mean a double bond; three would mean a triple bond. Letters such as N, O, Cl
and Br name atoms that must be shown. A `+` or `−` beside an atom shows its formal charge.

In the app, the formula and molecular weight appear in the status bar as you draw. If the formula
is unexpected, check each line end, corner, label and bond order. **Structure → Check Structure**
can point out valence errors and overlapping atoms, but it cannot know which molecule you meant.

A wedge or hashed bond shows how a bond points out of or behind the page. For an ordinary skeletal
structure, Penzene needs that information before it can assign an **R** or **S** stereocentre.
**Format → Stereo Labels** then displays the calculated labels. A flat drawing with no specified
stereochemistry has no R/S label to calculate. See [Show calculated stereo labels](stereochemistry.md)
for a worked example.

Now [draw ethanol yourself](first-molecule.md). Keep the editable `.penz` file when you finish;
export SVG or PDF when you need a picture for notes or a report.
