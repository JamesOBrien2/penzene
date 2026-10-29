# Browser preview feasibility

Assessment for [#713](https://github.com/JamesOBrien2/penzene/issues/713), 29 October 2026.
Penzene currently ships desktop applications and a Python module. There is no supported
browser build or hosted editor. This is a source and dependency assessment, not a successful
cross-compilation: Emscripten and a Qt WebAssembly kit were not available for this assessment.

## Finding

A browser preview using the existing Qt canvas is plausible, but requires a port rather
than another release packaging job. The first experiment should link Penzene's C++ chemistry
dependencies with Qt WebAssembly, before adapting the full application. Success of RDKit.js
alone would not establish that the editor can link or that its chemistry behaves correctly.

[Qt's WebAssembly documentation](https://doc.qt.io/qt-6.10/wasm.html) lists Widgets, GUI,
Network and SVG support. Qt requires an Emscripten version matched to the chosen Qt minor
release; its WebAssembly libraries are static. Nested event loops need Asyncify/JSPI or
asynchronous replacements. Files normally live in temporary memory, with separate APIs for
uploading and downloading local files. Clipboard access is permission-dependent and cannot
be assumed to carry Penzene's custom formats. Threaded builds require cross-origin isolation;
a single-threaded first experiment avoids that hosting requirement. Fonts need explicit
attention because a browser cannot assume the desktop's installed fonts.

## Dependency and source audit

| Area | Current implementation | Consequence for an experiment |
|---|---|---|
| Qt | `CMakeLists.txt` links GUI/SVG in the core and Widgets/Network/PrintSupport in the app. | Use a matching host Qt for build tools and a separate WebAssembly Qt kit for target libraries. Remove desktop printing from the initial preview; do not assume PrintSupport parity. |
| Chemistry | `penzene_core` directly links SmilesParse, FileParsers, Depictor, GraphMol, ChemReactions, Descriptors, RDInchiLib, CIPLabeler and DistGeomHelpers. | Cross-compile this dependency closure, including InChI and transitive dependencies, with compatible exception and thread settings. Native pixi libraries cannot be linked into Wasm. |
| ChemDraw files | `penzene_chemdraw` is built from vendored sources with Expat and Boost. The non-Apple/non-Windows CMake branch adds Boost.Locale and selects Linux Unicode conversion. | Audit this branch for Emscripten and cross-compile its dependencies. Dropping binary CDX would require making it optional throughout its callers, not just omitting one link library. |
| Build tools and data | CMake uses `lrelease`, `rcc` and a generated NMR resource; the Python module uses nanobind. | Keep Qt generators on the host, package required resources, and start with `PENZENE_PYTHON=OFF` and `PENZENE_TESTS=OFF`. Existing desktop tests cannot serve as browser runtime tests. |

[RDKit MinimalLib's README](https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/README.md)
documents its WebAssembly distribution and build script. Its
[CMake target](https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/CMakeLists.txt)
provides useful upstream evidence for compiling chemistry libraries to Wasm, including optional
InChI and reaction support. It exports a JavaScript-facing subset, however; Penzene's
`src/Chem.cpp` calls the C++ API directly. Replacing that bridge with RDKit.js would be a separate
architecture change, not a shortcut for preserving this editor.

The existing `cmake/build-rdkit.sh` is a native release build. Its choices about optional
RDKit features are a starting inventory, not a tested Emscripten recipe. The current
`pixi.toml` targets Linux, macOS and Windows only, and `CMakeLists.txt` has no WebAssembly branch.

## Browser behavior to adapt

| Desktop behavior | Initial preview boundary |
|---|---|
| Open, Save As and export use filesystem paths in `src/MainWindow.cpp`. | Adapt to user-selected uploads and explicit downloads. `.penz` must round-trip before advertising saving. |
| Autosave uses `QStandardPaths` and `QSaveFile`; settings and user styles also use local storage. | Treat the session as temporary until persistence, synchronization and reload recovery have been tested. Do not label an in-memory write as a durable save. |
| Dialogs and menus use `exec()`; `pubchem::fetch` in `src/Online.cpp` runs a nested `QEventLoop`. | Audit these entry points for asynchronous operation or a compatible Asyncify/JSPI build. Initially omit online lookup and update checks; later test CORS and failure handling. |
| Copy/paste offers chemical MIME data and platform-specific Office integration. | Prove text SMILES and internal selection copy/paste first. Offer file download for full-fidelity transfer; make no Office/OLE or custom clipboard round-trip promise. |
| Printing and export share the drawing engine. | Test SVG/PNG download first; validate PDF separately from native printing. |
| Labels use system fonts and the app has extensive keyboard/accessibility behavior. | Bundle a redistributable font, verify chemical glyphs, and test focus, shortcuts, text entry and accessibility in each target browser. |

## Smallest useful build spike

1. Pin one Qt minor, its matching Emscripten toolchain, and the repository's RDKit version.
   Build a tiny Wasm link test that parses SMILES, generates 2D coordinates, computes a formula
   and round-trips a MOL block. Exercise InChI, stereo and invalid-input exceptions explicitly.
2. Cross-compile the remaining core dependencies and link `penzene_core`. Render one document
   containing labels, charges, radicals and arrows. Record missing features rather than
   silently substituting another chemistry or drawing engine.
3. Bring up the existing canvas with sample data and undo/redo. Add upload/download adapters
   and audit modal operations before exposing the full menus. Keep desktop builds intact.
4. In current Chrome, Firefox and Safari, draw/edit a molecule, use the keyboard, round-trip
   `.penz`, and download SVG/PNG. Test denied clipboard access, reload data loss/recovery,
   stereo preservation and malformed uploads. Record download size, startup time and memory
   on named hardware; no size or performance estimate has been measured here.

Only a passing spike would justify a supported preview build and CI target. Hosting,
offline caching, long-term browser storage and desktop feature parity remain separate work;
this assessment adds no public service or new runtime dependency.
