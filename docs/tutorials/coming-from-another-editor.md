# Coming from another editor

Bring in a ChemDraw XML drawing, make an edit with familiar hotkeys, choose a drawing style, and
save it for another editor. This example changes phenol into anisole, with the methoxy group shown
as **OMe**.

Example files: {download}`editor-migration.cdxml <files/editor-migration.cdxml>` (phenol) and
{download}`editor-migration-finished.penz <files/editor-migration-finished.penz>` (the result).
These examples were generated in Penzene and contain no third-party drawing assets.

## 1. Open the drawing

In your other editor, save a copy as **ChemDraw XML (.cdxml)**, then choose **File → Open…** in
Penzene. You can also drag the file onto its window. Try the example `.cdxml` first: its status bar
should show **C₆H₆O**, phenol.

CDXML brings across molecules, abbreviations, arrows, text, basic shapes and lone pairs. Check your
own drawing after import, especially features specific to its original editor. Keep the original
file while you work; **File → Save As… → Penzene document** gives you a separate `.penz` working copy.
See [File formats](../files.md) for the supported formats and their limits.

## 2. Transfer the shortcuts you use

Penzene's *hotspot* is the atom or bond under the pointer, marked in green. It stays active when
you move the pointer away. Keys act on that hotspot; press Esc to clear it before choosing a tool
with a key. Case matters: lowercase `o` means oxygen, uppercase `O` means OMe.

The following ChemDraw-style hotkeys work in Penzene. Other editors may use different defaults;
[Keyboard shortcuts](../keys.md) is the complete Penzene list.

| Context | Keys | Penzene action |
|---|---|---|
| Atom | `c`, `n`, `o`, `l`, `b` | carbon, nitrogen, oxygen, chlorine, bromine |
| Atom | `O`, `y`, `K` | OMe, Boc, t-Bu |
| Atom | `1`, `2`, `a` | sprout a single bond, a carbonyl-containing group, or phenyl |
| Atom | `Enter` or `=` | type an element, abbreviation or label |
| Bond | `1`, `2`, `3` | single, double, triple; `2` again swaps a double bond's side |
| Bond | `w`, `h` | wedge and hash |
| Bond | `4`–`8` | fuse a ring of that size |
| No hotspot | `x`, `X`, `j`, `e`, `t`, Space | single bond, chain, benzene, arrow, text, select |

On an atom, `2` depends on its neighbours: a ring carbon gets C=O, while a terminal carbon can
sprout an acetyl group. Point carefully and use Undo if the result is different from what you meant.
In Penzene, `'` sets an atom-map number; **Format → Atom Numbers** shows ordinary atom indices.

For the example, point at phenol's **OH** until the oxygen lights up, then type uppercase `O`.
The label becomes **OMe** and the formula becomes **C₇H₈O**. **Structure → Expand Abbreviations**
draws out the methyl if you prefer it explicit; leave it abbreviated for the figure below.

## 3. Match the drawing style

Choose **Format → Drawing Style → RSC** for this example. **ACS 1996** and **JDP** are also
available. This sets the physical bond length, line widths, double-bond spacing and font; exports
use that size. An imported drawing initially keeps its label-to-bond proportions. Choosing a
Penzene style makes it consistent with your other Penzene figures.

```{image} ../_static/tutorials/editor-migration.svg
:alt: Anisole in the RSC drawing style with an OMe abbreviation
:width: 240px
```

**Edit → Preferences** sets the style for new documents. It does not restyle the document you
already have open. See [Drawing styles](../figures.md#drawing-styles) for the dimensions.

## 4. Save for further editing

Use **File → Save As…**, choose **ChemDraw XML (.cdxml)**, and save `anisole-return.cdxml`. Open
that copy in your other editor and check the chemistry, labels and layout. Saving the Penzene
working copy as `.penz` keeps the app's full drawing information, including multiple pages.

To send a figure instead, use **File → Export…** for SVG, PDF or PNG. **Edit → Copy** carries
the picture and chemistry formats; **Edit → Copy As → SMILES** gives structure text. Paste accepts
SMILES text as well as supported drawing formats. These commands are separate, so choose the one
your destination needs.

Continue with [Caffeine from the keyboard](keyboard.md) or [A publication figure](publication-figure.md).
