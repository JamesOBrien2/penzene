# A reaction scheme

Draw the synthesis of aspirin from salicylic acid and acetic anhydride, label the arrow, and export
it for a slide or a paper.

```{image} ../_static/tutorials/aspirin-synthesis.svg
:alt: Salicylic acid plus acetic anhydride, with H2SO4 over the arrow and 90 °C under it, gives aspirin
```

Example file: {download}`aspirin-synthesis.penz <files/aspirin-synthesis.penz>` (the finished scheme).

## 1. Put the molecules in

The quickest way is to paste reaction SMILES. Copy this line and choose **Edit → Paste** (⌘V or Ctrl+V)
in an empty document:

```text
OC(=O)c1ccccc1O.CC(=O)OC(C)=O>>CC(=O)Oc1ccccc1C(=O)O
```

Penzene lays out the reactants with a **+** between them, a reaction arrow, and the product. To
draw them yourself instead, draw each molecule ([Caffeine from the keyboard](keyboard.md) shows the
fast way), then choose a reaction arrow from **Arrows** and drag it between them.

## 2. Line it up

Choose **Arrange → Arrange Scheme**. The molecules sit on one line, centred on the arrow, evenly
spaced. Run it again after any change.

## 3. Reagents and conditions

Press Esc, then `t` (with an atom under the pointer, `t` types a label instead), or choose the text
tool from **Shapes**. Click just above the arrow and type `H2SO4`.
Digits after a letter are set as subscripts, so it shows as H₂SO₄. Click below the arrow and type
`90 °C`. **Edit → Copy As → Reaction SMILES** gives the scheme back as text, with salicylic acid and
acetic anhydride as reactants and aspirin as the product.

## 4. Export

- **File → Export…** and pick SVG, PNG or PDF. Exports come out at the drawing style's real size,
  black on clear or white (**Edit → Preferences** sets the background, PNG resolution and scale).
- **Edit → Copy** puts the scheme on the clipboard as vector PDF, PNG and SVG at once; paste it into
  Word, PowerPoint or Keynote.
- An exported SVG, PNG or PDF carries the drawing inside it. Open it in Penzene to edit it again.

From a script, the example file gives the image above:

```sh
penzene --render aspirin-synthesis.penz aspirin-synthesis.svg
```

More on schemes, arrows and pages: [Figures and exports](../figures.md).
