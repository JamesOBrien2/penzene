# Journal drawing styles

Choose a preset in **Format → Drawing Style**, Preferences, `--drawing-style`,
or Python's `Document.style`. **Edit Styles… → Duplicate** makes any preset
editable. All dimensions below are physical points at 100% export scale;
1 cm = 72 / 2.54 points. Canvas zoom and the editor preview do not alter them.

| Preset | Labels / captions | Bond | Line | Bold | Wedge | Margin | Hash | Double spacing | Chain / snap |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| ACS 1996 | Arial 10 / 10, normal | 14.4 | 0.6 | 2 | 3 | 1.6 | 2.5 | 18% | 120° / 30° |
| JDP | IBM Plex Sans Light 10 / normal 10 | 14.17 | 0.879 | 1.814 | 2.721 | 1.162 | 1.814 | 18% | 120° / 30° |
| RSC | Helvetica 7 / 7, normal | 12.2 | 0.5 | 1.6 | 1.6 | 1.25 | 1.8 | 20% | 120° / 15° |
| Nature | Helvetica 6 / 6, normal | 10.8 | 0.595276 | 1.559055 | 2.338583 | 1.190551 | 1.700787 | 18% | 120° / 30° |
| Angewandte | Arial 12 / 12, normal | 17 | 0.75 | 2.601563 | 3.902344 | 2 | 2.601563 | 18% | 120° / 30° |
| Angewandte (8 pt final) | Arial 8 / 8, normal | 11.333333 | 0.5 | 1.734375 | 2.601563 | 1.333333 | 1.734375 | 18% | 120° / 30° |
| SYNLETT/SYNTHESIS | Arial 6 / 6, normal | 10.119685 | 0.566929 | 1.388976 | 2.083465 | 0.907087 | 1.757480 | 18% | 120° / 30° |

The renderer keeps fractional font sizes. Bond spacing can be relative to the
actual untrimmed bond or an absolute point distance. Label margins are measured
against glyph outlines; the minimum centre clearance is zero in these presets,
so a guessed label radius does not override a publisher's margin.

## Verified sources and choices

Checked against the following sources on 16 October 2026:

- **ACS 1996:** Revvity's [published ACS document settings](https://support.revvitysignals.com/hc/en-us/articles/4408234173332-ChemDraw-When-I-print-the-structures-the-line-thickness-for-chemical-bonds-has-increased-considerably-How-do-I-get-it-back-to-drawing-finer-lines-for-the-chemical-bonds)
  and the supplied ACS Document 1996 stationery. This is a named stationery,
  not a claim that every ACS journal requires the same figure format.
- **JDP:** the existing JDPReport.cds values recorded with the original preset.
  This is a personal report style, not a publisher specification.
- **RSC:** the current [RSC author instructions](https://www.rsc.org/publishing/publish-with-us/publish-a-journal-article/crystengcomm)
  specify bond 12.2 pt, line 0.5 pt, bold and stereo width 1.6 pt, hash 1.8 pt,
  spacing 20%, Helvetica/Arial 7 pt, chain angle 120° and fixed angle 15°.
  The point values take precedence over their approximate centimetre equivalents.
  Margin 1.25 pt comes from the RSC stationery; the current text does not specify it.
- **Nature:** [chemical structures guide](https://www.nature.com/documents/nr-chemical-structures-guide.pdf),
  pages 1–5. The measured settings on page 4 specify bond 0.381 cm, line 0.021 cm,
  bold 0.055 cm, margin 0.042 cm, hash 0.06 cm, spacing 18%, chain 120°,
  and normal Arial/Helvetica 6 pt. Captions may be 6 or 6.5 pt; the default is 6 pt.
  Applying the preset chooses plain charges, places the charge before a right-hand
  H, raises R-group numbers, keeps implicit hydrogens visible and uses skeletal
  carbon labels and discrete aromatic bonds. These remain editable choices.
- **Angewandte:** the [current graphics guidelines](https://onlinelibrary.wiley.com/page/journal/15213773/homepage/graphics-guidelines)
  link to the [official ChemDraw template](https://onlinelibrary.wiley.com/pb-assets/assets/vch/angew-cds-1553091522227.zip).
  Its `angew.cds` contains bond 17 pt, line 0.75 pt, bold/hash 2.6015625 pt,
  margin 2 pt, spacing 18%, chain 120° and normal Arial 12 pt. **Angewandte**
  reproduces those original settings. The guide recommends approximately 7–8 pt
  at publication size: **Angewandte (8 pt final)** scales every physical template
  measurement by 8/12. Use that preset at 100%, or size the original template
  appropriately in the manuscript; do not reduce both.
- **SYNLETT/SYNTHESIS:** [SYNLETT's 2026 instructions](https://cdn.sanity.io/files/8egl0zzn/production/b1e2fddfba11be24d9129ef9b03b41179c1d6db9.pdf)
  and [SYNTHESIS's 2026 instructions](https://cdn.sanity.io/files/8egl0zzn/production/ca4bc9c5ca1a310515b6495a7e3fc760a3a5a04d.pdf),
  page 1, specify final-size bond 0.357 cm, line 0.020 cm, bold 0.049 cm,
  margin 0.032 cm, hash 0.062 cm, spacing 18%, chain 120° and Arial 6 pt.
  The older installed stationery requires a 69% reduction; Penzene uses the
  final-size table already. Keep export scale at 100% to avoid reducing twice.

Where the publisher does not specify a wedge width, Penzene follows ChemDraw's
[documented geometry](https://chemapps.stolaf.edu/iupac/cdx/sdk/properties/BoldWidth.htm):
the wide end is 1.5 times bold width. RSC explicitly specifies its own 1.6 pt
stereo width. Where a source is silent about snap angles or script placement,
the existing editable application defaults apply; they are not publisher mandates.

## Complete figure review

The presets implement the document measurements and supported label conventions.
Some requirements depend on the scientific content or the manuscript and need
review of the finished figure:

| Source requirement | How to meet it |
|---|---|
| Nature: consistent bond lengths, angles and weights; discrete aromatic bonds except suitable delocalized ligands | Apply the preset, use Clean where appropriate, and inspect distorted or perspective structures. Do not use bold bonds merely for emphasis. |
| Nature: heteroatoms labelled, consistent methyl/functional-group abbreviations, non-standard abbreviations defined | Use atom labels and expanded structures as appropriate. The preset does not rewrite molecular connectivity or the author's explicit labels. |
| Nature: dissociated counterions, clear radical/charge marks, consistent metal-complex notation | Use actual charges/radicals and disconnected counterions. Check organometallic examples and avoid unnecessary dative arrows. |
| Nature: perspectival aromatic ligands may need absolute double-bond spacing | Duplicate the preset and select **Absolute (pt)** spacing. Choose a distance that keeps the compressed ring legible. |
| Nature: stereochemistry indicated with solid/hashed wedges; omit unnecessary explicit H; preserve configuration on reorientation | Use stereochemical bond tools and Check Structure; inspect the configuration and necessary stereohydrogens. |
| Nature: reaction, equilibrium, resonance, retrosynthesis and electron-flow arrows used for their stated meaning | Choose the appropriate arrow tool and anchor electron flow to the intended bond/atom. |
| Nature: bold only for compound numbers/descriptors; italics follow chemical conventions | Use compound labels and the Format toolbar for appropriate text. Keep ordinary labels/captions normal. |
| Nature: exceptionally large structures may need shorter bonds while retaining line widths, and labels no smaller than 5 pt | Duplicate and change **Bond length** and **Label size** independently. Ordinary export scaling changes every measurement, so it is not the same operation. |
| Angewandte: final text approximately 7–8 pt and figure width suitable for the manuscript | Use the final-size preset or size the template deliberately; check the individual journal's 84/178 mm column layout and accepted formats. |
| SYNLETT/SYNTHESIS: 83/173 mm single/double column figures; captions in manuscript text, not inside the graphic | Choose the appropriate page/layout, export at final scale and place the figure caption in the manuscript. Reaction conditions and compound labels can remain in the drawing. |
| Source-specific file formats, resolution, colour and figure dimensions | Use the current journal instructions and [export guide](figures.md). These depend on the submitted figure type; a drawing preset cannot choose them universally. |

Styles do not embed licensed font files. Install the specified font; the editor
warns when it is unavailable. PDF and SVG keep the rendered lettering as outlines.
A substituted font cannot be made exact by exporting it afterwards. Use black
ink unless colour is needed and permitted by the journal. Review the exported
figure at its actual insertion size before submission.

## Customization and exchange

The editor provides separate bond dimensions, relative/absolute spacing, label
and caption fonts/sizes/weights/italics, margins, chain/snap angles, script
size/offsets, hydrogen/carbon choices, charge placement/circles and aromatic
circle defaults. Duplicate, save, import or export JSON to share these choices.
Saved drawings carry a snapshot, which remains editable without the local library.

PDF, SVG and PNG use the same renderer as the preview. ChemDraw files carry
supported document settings and caption typography; ChemDraw controls its own
wedge and script geometry. Use Penzene's vector exports to preserve custom
rendering exactly.

```fish
penzene --render aspirin.mol aspirin.svg --drawing-style Nature
penzene --render aspirin.mol aspirin.svg --drawing-style 'SYNLETT/SYNTHESIS'
```
