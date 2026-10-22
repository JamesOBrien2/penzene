# Drawing

## Tools and the hotspot

The rail on the left holds the tools in six groups: **Select** (select, rotate in 3D, eraser, colour), **Bonds** (single,
double, triple, wedge, hash, interaction, partial, chain), **Rings**, **Atoms** (element,
charges), **Arrows**, and **Shapes** (lines, boxes, ellipses, orbitals), plus a separate **Text** button.
Click **Text** (or press `t` with no hotspot) to add or edit text where you click. Return or
clicking away keeps it, Shift+Return starts a new line, and Esc drops the edit. Click a group to open its tools
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

```{image} _static/tutorials/caffeine-keys.gif
:alt: Caffeine drawn with one click and twelve hotkeys
:width: 420px
```

Caffeine from one click and twelve keys, step by step: [Caffeine from the keyboard](tutorials/keyboard.md).

## Atoms and labels

- Press **Enter**, `=` or `t` on an atom (or click it with the Text tool) to type a label: an element
  (`Br`, `NH2`), an abbreviation (`OMe`, `Boc`, `TBS`, `CO2Me`, …), a SMILES fragment (drawn out), or
  any text (`R1`, `X`, `MgEt`). Text that isn't chemistry is drawn as written and counts as a generic atom.
- Isotopes: type the mass number before the element (`13C`, `18OH`, `15NH2`), or `D` / `T` for
  deuterium and tritium (`d` on an atom makes it D). They're drawn ¹³C, carried through SMILES
  (`[13CH4]`), MOL and ChemDraw files, and counted in the masses.
- `x` and `r` on an atom label it X and R. `R1`, `R2`… export as MDL R-groups.
- **Stands for** (Atom Properties) says what a variable label means: `N, O, S` or `H, Me, OMe`. Penzene
  writes it under the molecule (X = N, O, S) and rewrites that line whenever a definition changes, so edit
  the definitions rather than the line. Elements only make the atom a query atom list in MOL files and
  Edit → Copy As → SMARTS; groups stay a note. Typing alternatives as the label (`N,O,S` or `N/O/S`)
  does both at once: the atom becomes X (or the next free letter) standing for them. For searches, the MDL codes `A`, `Q`, `M` (and `AH`, `QH`,
  `XH`, `MH`) can be typed as labels too.
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
heterocycles, macrocycles, ligands and reagents, and Haworth, Fischer and Newman projections. Haworth and Fischer drawings give the right stereo
when converted to SMILES. **Structure → Save Selection as Template** adds your own.

The **Lab** section has a TLC plate: a baseline, a dashed solvent front, and starting material,
co-spot and product lanes. Each spot shows its Rf beside it, worked out from where it sits between the
baseline and the front, and updated as you move a spot or the front (Ctrl/Cmd+click to move one part
of the plate on its own). Copy a spot to add one, drag its handles to draw a streak, and edit the
lane labels like any text.

```{image} _static/templates-light.png
:alt: The Templates panel open beside aspirin, showing the amino acids
:class: shot only-light
```

```{image} _static/templates-dark.png
:alt: The Templates panel open beside aspirin, showing the amino acids, dark theme
:class: shot only-dark
```

## Selecting and arranging

**Edit → Find Substructure…** (Ctrl/Cmd+F) selects every match of a SMILES or SMARTS pattern on the
current page. A selection supplies the initial pattern. Try `C=O` for carbonyls or `[#7,#8]` for
nitrogen and oxygen atoms. The status bar reports the match count; colour, delete, copy and arrange
commands act on the selected atoms and the bonds joining them. Matches inside abbreviations select
their visible label. Cancelling or entering an invalid pattern keeps the selection unchanged.

With the **Erase** tool, drag across atoms, bonds, arrows and text to mark them for deletion.
Release to erase everything crossed in one undo step, or press Escape to cancel. A click still
erases one object. Neighbours left without bonds are removed too, as with a click erase.

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
- **Structure → Brackets** puts square or round brackets around the selected atoms, with a subscript
  such as *n* and a superscript at the top right: ‡ for a transition state, or a charge such as •+ or 3−.
- **Structure → Variable Attachment** draws a bond across the edge of the selected atoms (a ring's,
  say) for a substituent on any one of them. Formula and mass count it once; MOL (V3000) and
  ChemDraw files keep the positions.

## Arrows, text, shapes and colour

### Saved drawing styles

Nature, Angewandte and SYNLETT/SYNTHESIS join the ACS 1996, JDP and RSC presets.
Their measurements and publisher sources are listed in [Journal drawing styles](journal-styles.md).

**Format → Circled Charges** surrounds single `+` and `−` charge signs with a circle,
including charges on skeletal carbons. Multiple charges such as `2+` and `2−`
use ordinary superscripts. Returning to a single charge restores its circle.
This is a saved, undoable document setting; it changes the drawing, not chemistry.
The default is off. Canvas, print and graphics exports use the same appearance.

**Format → Drawing Style → Edit Styles…** opens the style editor. Built-in styles
are read-only: **Duplicate** one, give it a unique name, adjust its measurements
in points, and **Save**. Bonds, Labels, Captions and Layout group the controls.
Set line, bold and wedge widths, hash spacing, relative or absolute multiple-bond
spacing, label margins and minimum clearance independently. Captions can use a
different font, size, weight and italic setting from atom labels. Chain angle,
fixed angle increment, script size/position, hydrogen visibility, charge position/circles,
R-group numbering and carbon labels are also adjustable.

The preview shows a metal complex, stereochemistry and a caption. **100%** and
**200%** preserve the physical scale across styles; **Fit** fills the preview only.
Changing the preview scale does not change exported measurements.
**Import** and **Export** exchange versioned JSON definitions; **Delete** removes
only a saved user style. Names are unique without regard to case.

Saved styles appear in the Drawing Style menu, Preferences, `--drawing-style`,
and Python's `drawing_styles()` and `Document.style`. All three use the same
`Penzene/styles` directory under the operating system's user data location.
Applying a custom style embeds a snapshot in the drawing, so changes to the
library do not change existing documents. Reapply the style to take a newer
definition. Font files are not embedded; install the named font on another
computer for the same lettering. The editor reports unavailable fonts. An embedded
style remains editable even when its library entry is absent or different; save
it under a new name to add it to the local library.

PDF, SVG and PNG use the same renderer as the preview. ChemDraw exchange carries
its supported document dimensions and caption typography. ChemDraw has its own
wedge, script and charge-placement rules, so use Penzene's vector exports when
those custom details must look identical.

The **Format toolbar** above the canvas sets the selection’s colour, shape or ring fill,
and text font, size in points, bold and italic. Mixed selections show **Mixed**. Each
change is one undo step. With nothing selected, choose the tool and set the format
for the next object drawn. Font controls apply to text; Fill applies to boxes, ellipses
and rings. **None** removes a fill. The Colour brush shares the toolbar’s colour;
**Automatic** restores the drawing’s ink. Separate shape fills and text formatting
are saved in `.penz` and in ChemDraw files.

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
  circle. Right-click a box or ellipse and choose **Filled** to fill it with its colour.
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

### Additional template sources

The heterocycles, macrocycles, ligands and reagents added to the library use chemical names and
molecular connectivity checked against the PubChem records below. Their drawings are generated
by RDKit within Penzene and distributed under Penzene’s GPL-3.0-or-later licence.

| Category | Template | PubChem |
|---|---|---|
| Heterocycles | Pyrazole | [CID 1048](https://pubchem.ncbi.nlm.nih.gov/compound/1048) |
| Heterocycles | Oxazole | [CID 9255](https://pubchem.ncbi.nlm.nih.gov/compound/9255) |
| Heterocycles | Isoxazole | [CID 9254](https://pubchem.ncbi.nlm.nih.gov/compound/9254) |
| Heterocycles | Thiazole | [CID 9256](https://pubchem.ncbi.nlm.nih.gov/compound/9256) |
| Heterocycles | Pyrazine | [CID 9261](https://pubchem.ncbi.nlm.nih.gov/compound/9261) |
| Heterocycles | Pyridazine | [CID 9259](https://pubchem.ncbi.nlm.nih.gov/compound/9259) |
| Heterocycles | 1H-1,2,4-triazole | [CID 9257](https://pubchem.ncbi.nlm.nih.gov/compound/9257) |
| Heterocycles | Benzofuran | [CID 9223](https://pubchem.ncbi.nlm.nih.gov/compound/9223) |
| Heterocycles | Benzothiazole | [CID 7222](https://pubchem.ncbi.nlm.nih.gov/compound/7222) |
| Heterocycles | Pyrrolidine | [CID 31268](https://pubchem.ncbi.nlm.nih.gov/compound/31268) |
| Macrocycles | 12-Crown-4 | [CID 9269](https://pubchem.ncbi.nlm.nih.gov/compound/9269) |
| Macrocycles | 15-Crown-5 | [CID 36336](https://pubchem.ncbi.nlm.nih.gov/compound/36336) |
| Macrocycles | Cyclen | [CID 64963](https://pubchem.ncbi.nlm.nih.gov/compound/64963) |
| Macrocycles | Cyclam | [CID 64964](https://pubchem.ncbi.nlm.nih.gov/compound/64964) |
| Ligands | 2,2'-Bipyridine | [CID 1474](https://pubchem.ncbi.nlm.nih.gov/compound/1474) |
| Ligands | 1,10-Phenanthroline | [CID 1318](https://pubchem.ncbi.nlm.nih.gov/compound/1318) |
| Ligands | Triphenylphosphine | [CID 11776](https://pubchem.ncbi.nlm.nih.gov/compound/11776) |
| Ligands | 1,2-Bis(diphenylphosphino)ethane (dppe) | [CID 74267](https://pubchem.ncbi.nlm.nih.gov/compound/74267) |
| Ligands | Ethylenediamine | [CID 3301](https://pubchem.ncbi.nlm.nih.gov/compound/3301) |
| Ligands | Acetylacetone | [CID 31261](https://pubchem.ncbi.nlm.nih.gov/compound/31261) |
| Reagents | 4-Dimethylaminopyridine (DMAP) | [CID 14284](https://pubchem.ncbi.nlm.nih.gov/compound/14284) |
| Reagents | Triethylamine | [CID 8471](https://pubchem.ncbi.nlm.nih.gov/compound/8471) |
| Reagents | N,N′-Dicyclohexylcarbodiimide (DCC) | [CID 10868](https://pubchem.ncbi.nlm.nih.gov/compound/10868) |
| Reagents | 1,1′-Carbonyldiimidazole (CDI) | [CID 68263](https://pubchem.ncbi.nlm.nih.gov/compound/68263) |
