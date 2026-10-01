# Drawing

## Tools and the hotspot

The rail on the left holds the tools in six groups: **Select** (select, rotate in 3D, eraser, colour), **Bonds** (single,
double, triple, wedge, hash, interaction, partial, chain), **Rings**, **Atoms** (element,
charges), **Arrows**, and **Shapes** (lines, boxes, ellipses, text). Click a group to open its tools
beside the rail. They stay open, so you can keep several on screen and drag each by its title or
any bare part of it, or drag a corner to resize it (the tools reflow to fit); close one with its **✕**. A new one opens
beside those already open, and Penzene reopens them where you left them; **View → Reset Tool Layout** puts them back. Clicking a group also
picks the tool you last used from it. Hover over a tool to see its key.

Point at an atom or bond and it becomes the **hotspot**, marked in green. It stays put when the mouse
moves off, so you can keep typing. Keys act on the hotspot: `1` sprouts a bond, `2`
a carbonyl, `a` a phenyl, `O` an OMe, and `+`/`−` change the charge. The arrow keys walk the
molecule atom → bond → atom. Every key is listed in [Keyboard shortcuts](keys.md) and under
**Help → Keyboard Shortcuts** (F1).

With no hotspot (press Esc), keys pick tools: `x` bond, `X` chain, `j` benzene, `e` arrow, `t` text,
Space select.

## Atoms and labels

- Press **Enter**, `=` or `t` on an atom (or click it with the Text tool) to type a label: an element
  (`Br`, `NH2`), an abbreviation (`OMe`, `Boc`, `TBS`, `CO2Me`, …), a SMILES fragment (drawn out), or
  any text (`R1`, `X`, `MgEt`). Text that isn't chemistry is drawn as written and counts as a generic atom.
- Isotopes: type the mass number before the element (`13C`, `18OH`, `15NH2`), or `D` / `T` for
  deuterium and tritium (`d` on an atom makes it D). They're drawn ¹³C, carried through SMILES
  (`[13CH4]`), MOL and ChemDraw files, and counted in the masses.
- `x` and `r` on an atom label it X and R. `R1`, `R2`… export as MDL R-groups.
- **Structure → Expand Abbreviations** draws them out in full.
- Atom Properties (`/` on an atom) sets the charge, map number, lone pairs, radical electrons and δ±.
- `:` cycles lone pairs, `*` toggles a radical dot, and the atom's context menu adds δ+ or δ−.
  Radicals are chemistry (one H fewer each); lone pairs and δ are only drawn.
- `.` adds an attachment point. `j` and `J` bond an η⁵-cyclopentadienyl or η⁶-benzene through the
  ring's centre.

## Bonds

- `1 2 3` set the order; `2` on a double bond moves its second line to the other side, and `l c r`
  place it left, centred or right.
- `w` / `h` wedge and hash (again to flip); `b` bold, `d` dashed, `y` wavy.
- `i` makes a bond an **interaction** (dotted: hydrogen bonds, contacts, coordination), and `p` / `P`
  a **partial** bond (dashed single or solid + dashed: bonds forming or breaking in a transition
  state). Neither counts as a covalent bond, so a drawn transition state keeps its reactants'
  formula.
- `f` brings a bond to the front: bonds it crosses are drawn with a gap.
- **Stereo groups** (enhanced stereo): right-click a stereocentre (or a selection) and pick
  **Stereo Group → Absolute**, **And n** or **Or n**. `&1` marks centres drawn as one of a mixture
  with their mirror image (racemic), `or1` centres that are one or the other, unknown which; centres
  sharing a number go together. The tag is drawn beside the centre and kept in MOL V3000 and ChemDraw
  files, and read from and written as CXSMILES (`C[C@H](N)C(=O)O |&1:1|`) by Copy as SMILES and `to_smiles()`;
  reaction SMILES leave them out.

## Rings and templates

Ring tools draw 3- to 8-membered rings and benzene: click empty space, an atom (spiro or attached),
or a bond (fused). On a bond, `a z v 4–8 9 0` fuse rings, `9` and `0` chair cyclohexanes.

**View → Templates** opens the template library: amino acids, sugars, nucleobases, common scaffolds,
and Haworth, Fischer and Newman projections. Haworth and Fischer drawings give the right stereo
when converted to SMILES. **Structure → Save Selection as Template** adds your own.

```{image} _static/templates-light.png
:alt: The Templates panel open beside aspirin, showing the amino acids
:class: shot only-light
```

```{image} _static/templates-dark.png
:alt: The Templates panel open beside aspirin, showing the amino acids, dark theme
:class: shot only-dark
```

## Selecting and arranging

- Drag to select, or Alt+drag from empty space to draw a loop round what you want; double-click selects a whole molecule. Drag a selection to move it, and drop an
  atom on another to merge them.
- The arrow keys nudge a selection (Shift for 10 points). Ctrl+arrow duplicates it across the next
  reaction arrow.
- **Rotate:** drag the round handle above a selection (Shift turns in 15° steps; Ctrl lands it square to the
  page or at 45°, however it started), Alt+drag,
  or Alt+←/→ in 15° steps. **Out of the page:** choose **Rotate in 3D** from the Select tool's
  flyout and drag a molecule, or use Shift+Alt+drag or Shift+Alt+arrows. Stereo is preserved.
- **Stretch and squash:** drag the handles around a selection (corners scale, edges stretch), or use
  **Arrange → Transform** for exact values.
- **Arrange → Align and Distribute**, **Center on Page** (for the selection or the whole
  drawing), **Flip**, and **Arrange Scheme** (lines up a reaction scheme,
  with agents centred over their arrows).
- **Arrange → Group** (Ctrl/Cmd+G) makes the selected molecules, arrows and text one object:
  clicking any part selects it all, and Arrange Scheme, Align and Flip move it as one. Ctrl/Cmd+click
  picks a single part inside a group. **Ungroup** (Shift+Ctrl/Cmd+G) splits it again. A selection
  shows one box per molecule (one per group); drag from anywhere inside a box to move it.
- **Structure → Brackets** puts square or round brackets, with a subscript such as *n*, around the
  selected atoms.
- **Structure → Variable Attachment** draws a bond across the edge of the selected atoms (a ring's,
  say) for a substituent on any one of them. Formula and mass count it once; MOL (V3000) and
  ChemDraw files keep the positions.

## Arrows, text, shapes and colour

- Arrows: reaction, no reaction (crossed), equilibrium, resonance, retrosynthesis, and curved (electron pair) and fishhook
  (single electron) arrows. Click a curved arrow again to flip its curve. While you draw a curved arrow, the atom
  or bond under each end lights up and the end settles on it; the arrow then moves with it when the structure is
  moved, turned or cleaned up.
- Select a single arrow to reshape it: drag either end, or drag the round handle at the top of a curve to make
  it deeper, shallower or bow the other way. Format → Arrowhead Size sets the heads of the selected arrows.
- Text: formulas get subscripts automatically (H2O → H₂O), and a charge at the end of a formula
  is set as a superscript (NH4+ → NH₄⁺, Cu2+ → Cu²⁺, [Fe(CN)6]3- → [Fe(CN)₆]³⁻). Where the digits
  could be either, mark the charge with ^: SO4^2- → SO₄²⁻.
- Lines (solid or dashed), boxes, rounded boxes and ellipses for grouping; Shift draws a square or
  circle.
- Orbitals, in **Shapes**: s, p, lobe and hybrid, each outlined, shaded or with a gradient. Click
  an atom to centre one on it (drag to point it), or click an orbital to restyle it; the colour
  tool colours it. New orbitals sit over the drawing: right-click for **Bring to Front**, **Bring
  Forward**, **Send Backward** and **Send to Back** (behind the bonds). ChemDraw files keep them,
  and their order.
- The colour tool paints atoms, bonds, arrows and text; the fill tool shades a ring.

## Pages, grid and rulers

- A drawing can have several pages, as tabs along the bottom of the window. **+** adds one; drag
  a tab to reorder it, double-click to rename it, and right-click (or the **Page** menu) to rename or delete it.
  Ctrl+PgDown and Ctrl+PgUp step through them. Each page has its own undo history.
- **Page → Move Selection To** moves the selection to another page. Copy and paste work between pages
  too.
- Pages are saved together in a `.penz` file. Other formats (MOL, ChemDraw) hold one page, so save
  as `.penz`, or export the page you want.
- **View → Grid** and **View → Rulers** measure the drawing at its final size: 5 mm grid squares,
  and rulers in centimetres from the page's corner. Both start on; turn either off there and
  Penzene remembers. On a trackpad, spread two fingers to zoom in and pinch to zoom out, and drag
  with two fingers to pan.
