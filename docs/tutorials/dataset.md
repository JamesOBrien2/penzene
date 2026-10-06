# A dataset of structures

Open a list of structures, look at each one's properties, and export a descriptor table: the
everyday side of a dataset, with what Penzene does today.

Example file: {download}`drugs.smi <files/drugs.smi>`, six drugs as SMILES, each followed by its name:

```text
CC(=O)Oc1ccccc1C(=O)O aspirin
Cn1c(=O)c2c(ncn2C)n(C)c1=O caffeine
CC(C)Cc1ccc(cc1)C(C)C(=O)O ibuprofen
CC(=O)Nc1ccc(O)cc1 paracetamol
CN1CCC[C@H]1c1cccnc1 nicotine
COc1ccc2cc(ccc2c1)[C@H](C)C(=O)O naproxen
```

## 1. Open it as a grid

**File → Open…** and choose `drugs.smi`. Each line becomes a structure, laid out in a grid. `.sdf`
and `.inchi` files open the same way.

```{image} ../_static/tutorials/drugs-grid.svg
:alt: Aspirin, caffeine, ibuprofen, paracetamol, nicotine and naproxen in a three-by-two grid
```

## 2. Look at one

Click a structure to select it and open **View → Properties Panel** (Ctrl+I). It shows the formula,
masses, elemental analysis, cLogP, TPSA, hydrogen-bond donors and acceptors, and the Lipinski and Veber
checks for the selection. Click another structure and the panel follows. Clean up a drawing,
recolour it, or correct it here like any other.

## 3. Export a descriptor table

**File → Export Descriptors…** writes a CSV with one row per molecule in the selection, or on the
page when nothing is selected. The command line writes the same table without opening a window:

```sh
penzene --descriptors drugs.smi --columns name,formula,mw,clogp,tpsa,lipinski_violations
```

```text
name,formula,mw,clogp,tpsa,lipinski_violations
aspirin,C9H8O4,180.16,1.31,63.60,0
caffeine,C8H10N4O2,194.19,-1.03,61.82,0
ibuprofen,C13H18O2,206.28,3.07,37.30,0
paracetamol,C8H9NO2,151.16,1.35,49.33,0
nicotine,C10H14N2,162.24,1.85,16.13,0
naproxen,C14H14O3,230.26,3.04,46.53,0
```

Leave out `--columns` for every column; add `--out drugs.csv` to write a file. The columns are listed
under [Descriptor tables](../cli.md#descriptor-tables). A line that can't be read keeps its row,
with the reason in the `error` column.

## 4. Save the set

**File → Save As…** with `.sdf` writes one record per molecule, for other software. `.penz` keeps
the grid as drawn.

To render each structure to its own image, see [Batch rendering from the command line](cli-batch.md).
