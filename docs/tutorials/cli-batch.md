# Batch rendering from the command line

`penzene --render` draws structures without opening a window. Here it turns a folder of MOL files, a
SMILES library and a few SMILES strings into images, for a slide deck, a website or a supporting
information file.

Example files: {download}`aspirin.mol <files/molecules/aspirin.mol>`,
{download}`caffeine.mol <files/molecules/caffeine.mol>`,
{download}`paracetamol.mol <files/molecules/paracetamol.mol>` (put them in a `molecules` folder),
and {download}`drugs.smi <files/drugs.smi>`.

On macOS the program is inside the app: use `/Applications/Penzene.app/Contents/MacOS/penzene`, or add
an alias such as `alias penzene=/Applications/Penzene.app/Contents/MacOS/penzene`.

## Every MOL file in a folder

A shell glob passes them all at once; `--out` names the folder and `--format` the kind of image:

```sh
penzene --render molecules/*.mol --out out --format png
```

```text
out/aspirin.png
out/caffeine.png
out/paracetamol.png
```

Each image is named after its input file. The paths written are printed one per line, ready to pipe
into another command.

## A library, and SMILES strings

A `.smi` or `.sdf` file gives one image per record, named after the record (a `.smi` line's second
column, an SDF title). SMILES strings can go on the same command line:

```sh
penzene --render drugs.smi "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --out out --drawing-style RSC --clean
```

```text
out/aspirin.svg
out/caffeine.svg
out/ibuprofen.svg
out/paracetamol.svg
out/nicotine.svg
out/naproxen.svg
out/structure.svg
```

- SVG is the default format.
- A SMILES string has no name, so it becomes `structure.svg`.
- `--drawing-style RSC` draws in the RSC style (`ACS 1996` is the default; `JDP` is the third).
- `--clean` lays every structure out afresh instead of using the file's coordinates.

## One file with everything

Give an output path instead of `--out`, and the records are laid out as a grid in a single image:

```sh
penzene --render drugs.smi drugs.pdf
```

```{image} ../_static/tutorials/drugs-grid.svg
:alt: The six structures of drugs.smi in a grid
```

Every option, and the descriptor tables: [Command line](../cli.md).
