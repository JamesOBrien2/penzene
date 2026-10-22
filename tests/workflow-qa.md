# Manual user workflows

Use these tasks to find bugs, missing capabilities and confusing instructions that unit tests
cannot catch. Give a participant the goal and the linked documentation, but let them choose the
controls. Record where they hesitate or need a hint. Run the same task yourself afterwards to
separate a reproducible defect from a documentation or discoverability problem.

Start with the latest stable download linked from the [installation page](../docs/install.md).
Record its version, OS, display scale, and whether this is a fresh install. Use a separate test
account or VM for fresh-preference tests; do not reset anyone's working settings. Keep the sample
files in their own folder so a Save As cannot overwrite them.

| ID | User and task | Start with | Completion check |
|---|---|---|---|
| W1 | New undergraduate: install, draw ethanol, save, reopen and export a figure | [Reading drawings](../docs/tutorials/reading-structures.md), [first molecule](../docs/tutorials/first-molecule.md) | They find the correct installer, make C₂H₆O / MW about 46.07, reopen an editable `.penz`, and obtain a sharp SVG without help. Repeat on macOS and Windows. |
| W2 | Student: identify a stereocentre and show its assignment | [Stereo labels](../docs/tutorials/stereochemistry.md), the two supplied glyceraldehyde files | **Format → Stereo Labels** shows R for the first drawing and S for its mirror; the state survives save/reopen and export. A flat unspecified centre is not given an invented R/S label. Note whether the participant finds the control without a hint. |
| W3 | Synthetic chemist: make and revise a reaction scheme | [A reaction scheme](../docs/tutorials/reaction-scheme.md) | Aspirin reaction SMILES imports; reactants, arrow, `H2SO4` and `90 °C` remain legible after Arrange Scheme, a later edit, save/reopen, and SVG/PDF export. |
| W4 | Student writing a report: move a figure into Word or PowerPoint | [Office and LaTeX](../docs/tutorials/office-latex.md), `caffeine-figure.penz` | SVG stays sharp and correctly proportioned at final size; the `.penz` source remains editable. On Windows with the installer, test the Penzene Drawing Object separately from a static picture. |
| W5 | ChemDraw user: continue a collaborator's drawing | [Coming from another editor](../docs/tutorials/coming-from-another-editor.md), `editor-migration.cdxml` | CDXML opens with labels and layout intact, accepts an edit, saves as `.penz` and CDXML, then reopens without losing the edit. If ChemDraw is available, check the saved CDXML there too. |
| W6 | Lab researcher: prepare supporting information | [Supporting information](../docs/tutorials/supporting-information.md), its sample files | Properties and descriptor CSV have the expected molecules and units; predicted NMR text is clearly marked as predicted. Compare the exported files to the tutorial's supplied examples. |
| W7 | Computational chemist: render a set of structures | [Dataset](../docs/tutorials/dataset.md), [Python batch](../docs/tutorials/python-batch.md) | A `.smi` library opens as a grid, descriptor CSV has one row per readable record, and the example script makes one figure per input. Record any skipped or malformed record visibly rather than silently accepting it. |
| W8 | Inorganic chemist: inspect metal complexes | `tests/data/metals.smi` (especially cisplatin and ferrocene), [Drawing](../docs/drawing.md#bonds) | Donor labels and hydrogens remain visible, formulas are correct, Clean keeps a readable geometry, and a save/reopen does not change coordination or charges. |

For W1–W5, try both a Mac and a Windows computer if available. Use at least two people who are
new to Penzene for W1 and W2, and one experienced chemistry-drawing user for W3–W6. W7 and W8
need someone familiar with the input chemistry. Time the first attempt, note every hint given,
and ask the participant what they expected before showing a workaround.

## Recording findings

Create one note per observation with the workflow ID, app version, OS, goal, steps taken,
expected result, actual result, whether it repeats, and a screenshot or minimal `.penz`/input
file where useful. Mark it **bug** (wrong or lost result), **missing capability** (task cannot be
done), **discoverability** (capability exists but could not be found), or **docs** (instructions
were unclear or wrong). Reproduce bugs with the smallest safe example before opening a GitHub
issue; avoid posting private lab structures. Link related observations to one root-cause issue
instead of creating duplicates.

Triage in this order: crashes/data loss or wrong chemistry; inability for a new user to finish
W1 or W2; broken common export or file exchange; then friction and polish. A missing R/S label
for a drawing with specified stereo is a chemistry bug. A participant failing to find the
existing menu toggle is a discoverability finding; it may justify a visible button if repeated.
Retest each accepted fix on its original workflow and on the other platform before calling it
done. Keep the resulting bug and feature list separate from the release notes until fixes ship.
