# Batch rendering from Python

The `penzene` module draws with the same engine as the app. Here a short script reads a CSV, draws
one figure per row and prints each formula; then the same drawings appear inline in Jupyter.

Example files: {download}`compounds.csv <files/compounds.csv>` and
{download}`render_csv.py <files/render_csv.py>`.

```text
name,smiles
aspirin,CC(=O)Oc1ccccc1C(=O)O
caffeine,Cn1c(=O)c2c(ncn2C)n(C)c1=O
ibuprofen,CC(C)Cc1ccc(cc1)C(C)C(=O)O
paracetamol,CC(=O)Nc1ccc(O)cc1
```

## The script

```{literalinclude} files/render_csv.py
:language: python
```

Run it next to `compounds.csv`:

```sh
python render_csv.py
```

```text
aspirin C9H8O4 180.16
caffeine C8H10N4O2 194.19
ibuprofen C13H18O2 206.28
paracetamol C8H9NO2 151.16
```

It writes `figures/aspirin.svg`, `figures/caffeine.svg` and so on. Here is `figures/ibuprofen.svg`:

```{image} ../_static/tutorials/python-ibuprofen.svg
:alt: Ibuprofen in the RSC style
```

- `pz.from_smiles` lays the structure out; `clean()` lays it out again with RDKit, which matters for
  structures read from files with poor coordinates.
- `doc.style` takes `"ACS 1996"`, `"RSC"` or `"JDP"`.
- `export` writes `.svg`, `.png` or `.pdf` by the file's extension, exactly as the app exports. A PNG
  comes out at 300 dpi; `doc.to_png(dpi=600)` returns the bytes at another resolution.

## From an SDF instead

`pz.read("library.sdf")` opens every record laid out as one grid, as the app does. To draw them
one at a time, read the records with RDKit and pass each SMILES to `pz.from_smiles`:

```python
from rdkit import Chem
import penzene as pz

for mol in Chem.SDMolSupplier("library.sdf"):
    if mol is None:
        continue                                 # a record RDKit can't read
    doc = pz.from_smiles(Chem.MolToSmiles(mol))
    doc.export(f"figures/{mol.GetProp('_Name')}.svg")
```

## In Jupyter

A `Document` displays itself: end a cell with one and the drawing appears below it, as SVG.

```python
import csv
import penzene as pz

with open("compounds.csv", newline="") as f:
    docs = {row["name"]: pz.from_smiles(row["smiles"]) for row in csv.DictReader(f)}
docs["caffeine"]
```

To show several from one cell, use `IPython.display.display(doc)` for each.

More: [Python](../python.md) and the [Python API reference](../python-api.md).
