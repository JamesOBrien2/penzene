# Supporting information and descriptor tables

Copy a molecule's calculated properties, export a table for a dataset, and use the NMR panel's
predicted SI line as a starting point for reporting measured data.

Example files:

- {download}`si-compounds.smi <files/si-compounds.smi>`: aspirin and ethanol.
- {download}`si-compounds.csv <files/si-compounds.csv>`: the GUI descriptor export for both molecules.
- {download}`si-aspirin-properties.txt <files/si-aspirin-properties.txt>`: copied Properties output.
- {download}`si-ethanol-predicted-nmr.txt <files/si-ethanol-predicted-nmr.txt>`: a prediction, not an experimental spectrum.

## 1. Select one complete molecule

Open `si-compounds.smi` with **File → Open…**. The two structures appear together. Choose **Select**
and double-click aspirin to select the entire molecule; selecting only an atom or part of it would
give properties for that fragment. Open **View → Properties Panel** (⌘I or Ctrl+I).

```{image} ../_static/tutorials/si-properties.png
:alt: Aspirin selected beside ethanol, with the Properties panel showing its formula, masses and descriptors
:width: 850px
```

The formula is **C₉H₈O₄**, average molecular weight **180.16**, and neutral monoisotopic mass
**180.0423**. The exact mass here is for the drawn molecule, rather than a protonated or sodiated
ion; the Mass Spec panel supplies those ion calculations.

## 2. Copy the properties

Click **Copy as Text** in Properties and paste into your notes or supporting-information draft.
The example's output is:

```{literalinclude} files/si-aspirin-properties.txt
:language: text
```

These are calculated values. cLogP is Wildman–Crippen logP, TPSA is in Å², and elemental percentages
come from the drawn formula. Keep those labels and units when you paste, and report experimental
elemental analysis or mass measurements separately. The Lipinski and Veber results are the checks
shown in [Chemistry](../chemistry.md), not measured biological activity.

## 3. Export every molecule as a CSV

Press Esc to clear the selection, then choose **File → Export Descriptors…** and save
`si-compounds.csv`. With nothing selected, each molecule on the current page gets a row. With a
selection, only its molecules or fragments are exported, so select complete structures deliberately.

The GUI uses generated row names such as `structure-1` and `structure-2` in the downloadable
example. Its `smiles`, `inchikey` and `formula` columns identify the chemistry. The `error` column
is empty for these valid structures; check it when importing a larger library.

For original record names and a chosen set of columns, run the command line beside the `.smi` file:

```sh
penzene --descriptors si-compounds.smi --columns name,formula,mw,exact_mass,clogp,tpsa --out si-selected-columns.csv
```

The resulting table is:

```text
name,formula,mw,exact_mass,clogp,tpsa
aspirin,C9H8O4,180.16,180.0423,1.31,63.60
ethanol,C2H6O,46.07,46.0419,-0.00,20.23
```

The small negative cLogP rounded to `-0.00` is effectively zero at this precision. Open the CSV in
your spreadsheet or analysis tool; preserve the identifiers when joining it to measurements.
[Descriptor tables](../cli.md#descriptor-tables) lists every column. Include `error` when choosing
columns for a dataset whose structures have not already been checked.

## 4. Copy an NMR line

Select the complete ethanol molecule and open **View → NMR Panel**. Choose **¹H** in the nucleus
selector. The panel follows the selected molecule; without a selection it uses the largest one.
Point at a predicted stick to light the atoms it belongs to, or focus the spectrum and press Left
and Right. A `~` indicates a weaker prediction match.

```{image} ../_static/tutorials/si-nmr.png
:alt: Ethanol selected with its predicted proton NMR sticks in the NMR panel
:width: 850px
```

Click **Copy SI Line**. The example produces:

```{literalinclude} files/si-ethanol-predicted-nmr.txt
:language: text
```

This is labelled **predicted**. It helps associate signals with the drawn groups and supplies a
format for an SI paragraph. Predictions do not account for solvent or stereochemistry, and can
miss distinctions such as diastereotopic protons. Replace the shifts, multiplicities and integrals
with your measured assignments, and add the actual solvent and spectrometer frequency. Keep the
prediction label until the line contains measured data.

More: [NMR prediction](../chemistry.md#nmr-prediction) and
[A dataset of structures](dataset.md).
