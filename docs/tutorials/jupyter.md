# Penzene documents in Jupyter

Use a notebook to show a drawing, edit an independent copy, and save figures and editable documents.
This example turns phenol into anisole and reopens the exported drawing. It goes beyond batch
rendering: the document stays available for edits between cells.

Download {download}`penzene-documents.ipynb <files/penzene-documents.ipynb>` and run it from top to
bottom. It creates its own input structures, so no data files are needed. The finished editable
example is {download}`notebook-anisole.penz <files/notebook-anisole.penz>`.

## 1. Use the notebook's Python environment

You need Python 3.12 or newer and Penzene. Install JupyterLab in that same environment:

```sh
python -m pip install penzene jupyterlab
python -m jupyter lab
```

Open the notebook and choose a plain Python kernel. If a notebook says it cannot import `penzene`,
install the package into that kernel's environment, then restart the kernel. Penzene includes Qt;
avoid a kernel that has already loaded PyQt/PySide or enabled a Qt GUI event loop. The example's
editing calls use the public API available in Penzene 1.9.

## 2. Let a document display itself

```python
from pathlib import Path
from IPython.display import display
import penzene as pz

original = pz.from_smiles("Oc1ccccc1")
original.style = "RSC"
original
```

Ending a cell with the `Document` displays an SVG drawing inline. It is the same rendering engine
as the desktop app, and does not need an extra image wrapper.

```{image} ../_static/tutorials/notebook-phenol.svg
:alt: Phenol rendered inline by a Penzene Document
:width: 240px
```

The document also exposes calculated properties:

```python
print(original.formula)
print(f"{original.mw:.2f}")
print(f"{original.exact_mass:.4f}")
```

For phenol these are `C6H6O`, `94.11` and `94.0419`.

## 3. Edit a copy

```python
edited = pz.from_json(original.to_json())
oxygen = next(i for i, atom in enumerate(edited.atoms) if atom.symbol == "O")
edited.hotkeys(oxygen, "O")
edited.clean()
edited
```

Uppercase `O` is the same OMe hotkey as in the desktop editor. Phenol becomes anisole, **C₇H₈O**;
the original remains **C₆H₆O**. Copying through `.penz` JSON keeps the drawing settings and gives
an independent document. `edited = original` would make two names for one mutable object.

```{image} ../_static/tutorials/notebook-anisole.svg
:alt: Anisole after replacing phenol's OH with OMe and cleaning the drawing
:width: 240px
```

Atom properties are read-only views; use the document's editing methods. Find indices from the
current document after adding or removing atoms. The example finds oxygen by its symbol instead
of assuming a fixed index. See [Python API reference](../python-api.md) for the available methods.

Only the last value in a cell displays automatically. To compare drawings in one cell:

```python
display(original)
display(edited)
```

## 4. Save and export

```python
output = Path("notebook-output")
output.mkdir(exist_ok=True)
edited.save(str(output / "anisole.penz"))
edited.export(str(output / "anisole.svg"))
edited.export(str(output / "anisole.pdf"))
(output / "anisole-600dpi.png").write_bytes(edited.to_png(dpi=600))
```

`.penz` keeps the editable source; SVG and PDF are vector figures, and PNG records the requested
resolution. The notebook writes these files into its own `notebook-output` folder. Rerunning the
cell replaces those example outputs.

## 5. Reopen and continue editing

```python
reopened = pz.read(str(output / "anisole.penz"))
assert reopened.to_smiles() == edited.to_smiles()
reopened
```

Penzene's exported SVG also carries its editable drawing, so `pz.read` can read it back. You can
open `anisole.penz` in the desktop app, edit and save it, then call `pz.read` again in the notebook.
An object already loaded in Python does not refresh automatically when the file changes.

The notebook ends with an exercise: replace phenol's OH with F in another independent copy,
predict the formula, and check the answer in the following cell. Its full example was executed
from a fresh kernel, including the inline displays, independent edits and export/reopen checks.

More: [Python](../python.md) and [Batch rendering from Python](python-batch.md).
