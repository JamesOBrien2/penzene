# Journal drawing styles

Choose a preset under **Format → Drawing Style**, in Preferences, with
`--drawing-style`, or through Python's `Document.style`. New presets are read-only;
duplicate one in **Edit Styles…** to adapt it. Measurements below are physical
points, before any additional export scaling; 1 cm = 72 / 2.54 points.

| Preset | Font | Label | Bond | Line | Bold | Margin | Hash | Double-bond spacing |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Nature | Helvetica | 6 | 10.8 | 0.595276 | 1.559055 | 1.190551 | 1.700787 | 18% |
| Angewandte | Arial | 12 | 17 | 0.75 | 2.601563 | 2 | 2.601563 | 18% |
| SYNLETT/SYNTHESIS | Arial | 6 | 10.119685 | 0.566929 | 1.388976 | 0.907087 | 1.757480 | 18% |

## Sources

Verified on 13 October 2026:

- [Nature's chemical structures guide](https://www.nature.com/documents/nr-chemical-structures-guide.pdf),
  page 4: bond 0.381 cm, line 0.021 cm, bold 0.055 cm, margin 0.042 cm,
  hash 0.06 cm, spacing 18%, and Arial or Helvetica 6 pt. This preset uses Helvetica.
- [Angewandte's current graphics guidelines](https://onlinelibrary.wiley.com/page/journal/15213773/homepage/graphics-guidelines)
  link to the [official ChemDraw template](https://onlinelibrary.wiley.com/pb-assets/assets/vch/angew-cds-1553091522227.zip).
  Its `angew.cds` contains bond 17 pt, line 0.75 pt, bold/hash 2.6015625 pt,
  margin 2 pt, spacing 18%, and Arial 12 pt. The preset follows that template.
  The current guidelines also recommend approximately 7–8 pt text at printed
  size. Check the final figure size; the historical 60% production reduction is
  not applied automatically by this preset.
- [SYNLETT's 2026 instructions](https://cdn.sanity.io/files/8egl0zzn/production/b1e2fddfba11be24d9129ef9b03b41179c1d6db9.pdf)
  and [SYNTHESIS's 2026 instructions](https://cdn.sanity.io/files/8egl0zzn/production/ca4bc9c5ca1a310515b6495a7e3fc760a3a5a04d.pdf),
  page 1: the final-size settings are bond 0.357 cm, line 0.02 cm, bold 0.049 cm,
  margin 0.032 cm, hash 0.062 cm, spacing 18%, and Arial 6 pt. Their older
  preinstalled template requires 69% reduction; this preset already uses the
  final-size table, so keep export scale at 100% to avoid shrinking twice.

The publishers do not specify a separate wedge width in those settings. Penzene
uses the existing ACS wedge-to-bold ratio of 2.25. The minimum label clearance
uses the existing nominal half-glyph width plus the stated margin. Lettering
depends on installed fonts; styles do not include font files.

## Presentation checks

These presets set measurements and typography. They do not enforce every journal
rule, change chemistry, choose a page width, or rewrite a figure. Review the
publisher's current requirements and the exported figure before submission.
Nature's guide discourages circled charges; leave **Circled Charges** off for
Nature figures. Applying its style preserves your chosen document setting.

```sh
penzene --render aspirin.mol aspirin.svg --drawing-style Nature
penzene --render aspirin.mol aspirin.svg --drawing-style 'SYNLETT/SYNTHESIS'
```
