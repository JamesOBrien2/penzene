# Chemistry

Penzene keeps the drawing and hands its chemistry to RDKit. Everything below updates as you draw.

```{image} _static/properties-light.png
:alt: The Properties panel for aspirin: formula, masses, cLogP, TPSA and the Lipinski and Veber checks
:class: shot only-light
```

```{image} _static/properties-dark.png
:alt: The Properties panel for aspirin: formula, masses, cLogP, TPSA and the Lipinski and Veber checks, dark theme
:class: shot only-dark
```

::::{container} features-grid

```{feature} scale
:title: Formula and mass
The status bar shows the formula, molecular weight and exact mass of the selection, or of the whole
drawing. **View → Properties Panel** (Ctrl+I) adds elemental analysis, cLogP, TPSA, H-bond donors
and acceptors, rotatable bonds, heavy atoms, and Lipinski and Veber checks, with **Copy as Text**
for a supporting-information table. **File → Export Descriptors…** writes them as a CSV, one row
per molecule ([descriptor tables](cli.md#descriptor-tables)).
```

```{feature} chart-bar
:title: Isotope patterns
**View → Mass Spec Panel** draws the isotope pattern of the selection, or of the whole drawing, as
a stick spectrum for [M], [M+H]⁺, [M+Na]⁺ or [M−H]⁻, from natural isotope abundances. The main
peaks are labelled with their m/z (aspirin [M+H]⁺ at 181.0495, with its ¹³C peak at 182.0529). A
drawn isotope such as ¹³C counts as that isotope only.
```

```{feature} atom
:title: Implicit hydrogens and valence
Labels get their hydrogens (OH, NH₂). An atom with too many bonds is drawn in red.
```

```{feature} rotate-3d
:title: Stereochemistry
Wedges and hashes set stereocentres; E/Z comes from the drawing. **Format → Stereo Labels**
shows CIP (R)/(S) and (E)/(Z).
```

```{feature} checklist
:title: Check Structure
The Structure menu lists valence errors, unknown labels, overlapping atoms, stereocentres without a
wedge, and wedges on atoms that aren't stereocentres. Click one to select it. Turn on **Check structures
before export and copy** in Preferences to see the list, with Export Anyway or Cancel, before a figure goes out.
```

```{feature} world-search
:title: Look Up on PubChem
**Structure → Look Up on PubChem** opens the selection's (or the drawing's) PubChem page in your browser,
matched by InChIKey. Penzene itself never goes online.
```

```{feature} wand
:title: Clean
Ctrl+Shift+K lays a structure out afresh, keeping its stereo and bond styles.
```

```{feature} letter-h
:title: Hydrogens
**Structure → Add or Remove Explicit Hydrogens**. **Format → Carbon Labels** changes how carbons are shown.
```

```{feature} hexagon
:title: Aromatic circles
The Format menu draws benzene-like rings with a circle, for the whole drawing or for selected rings.
```

```{feature} hash
:title: Numbers
**Format → Atom Numbers** shows atom indices; `'` on an atom sets its reaction atom-map number, which
goes into SMILES and MOL.
```

```{feature} atom
:title: Predicted NMR shifts
**Format → Predicted NMR Shifts** writes each carbon's predicted ¹³C shift beside it, and the ¹H shift
of any atom with hydrogens in brackets, in ppm. A `~` marks a weaker match. **View → NMR Panel** draws
the selection's (or the drawing's) predicted ¹³C or ¹H spectrum, one stick per set of equivalent atoms,
¹H sticks split first order (s, d, t, q, m) by the H on neighbouring carbons; point at a stick, or press
Left and Right, to light its atoms, or point at an atom to find its stick. **Copy SI Line** copies the prediction as a supporting-information
line ("1H NMR (predicted) δ 3.69 (q, 2H), 1.22 (t, 3H)."), to replace with measured values. See [NMR prediction](#nmr-prediction).
```

```{feature} world-search
:title: Names (online)
**File → Import → Name** turns a name into a structure, and **Edit → Copy As → IUPAC
Name** copies a structure's IUPAC name. Both look the compound up on PubChem, so they need an
internet connection and only work for compounds PubChem knows.
```

```{feature} copy
:title: Copy as
SMILES, InChI, InChIKey and reaction SMILES (Edit → Copy As).
```

::::

:::{note}
**What isn't chemistry.** Interaction and partial bonds, lone pairs, δ, brackets, arrows, shapes and
text are drawn only. Free-text labels are generic atoms (`*` in SMILES).
:::

## NMR prediction

Shifts are looked up by HOSE code in a table built from [nmrshiftdb2](https://nmrshiftdb.nmr.uni-koeln.de),
an open database of assigned spectra. A HOSE code describes an atom's surroundings sphere by sphere
(Bremser's scheme, written as the CDK writes it). Penzene uses up to four spheres. When the table has
no match at four, it tries three, then two, then one. The shift given is the median of every shift
recorded for that code. A `~` marks a match of two spheres or fewer.

¹H shifts are averaged over an atom's hydrogens, so diastereotopic CH₂ protons show as one value.
Solvent and stereochemistry aren't taken into account.

On one compound in ten held out of the table (by InChIKey), the mean absolute errors were:

| Spheres matched | ¹³C (ppm) | ¹H (ppm) | Share of ¹³C shifts |
|---|---|---|---|
| 4 | 1.1 | 0.15 | 54% |
| 3 | 2.5 | 0.31 | 24% |
| 2 | 4.9 | 0.53 | 19% |
| 1 | 10.8 | 0.94 | 3% |
| All | 2.4 | 0.27 | |

:::{note}
**Data licence.** Contains information from nmrshiftdb2 (www.nmrshiftdb.org), which is made available
here under the [nmrshiftdb2 Database License](https://nmrshiftdb.nmr.uni-koeln.de/nmrshiftdbhtml/nmrshiftdb2datalicense.txt).
The table (`resources/nmr/hose.tar.xz`) is a derivative database under that licence, not Penzene's
GPL-3, and a copy of the licence ships beside it. `cmake/nmr-table.py` rebuilds the table from
nmrshiftdb2's public export using Penzene's own HOSE codes. Drawings and exports that show predicted
shifts carry this notice under the structure.
:::
