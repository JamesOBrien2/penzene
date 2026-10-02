# Caffeine from the keyboard

Penzene's hotkeys put a structure together without changing tools. Point at an atom or bond to make
it the *hotspot* (it lights up), then type a key. Here is caffeine with one click and twelve keys.

```{image} ../_static/tutorials/caffeine-keys.gif
:alt: Caffeine drawn step by step: a cyclohexane, a fused five-membered ring, two C=O groups, four nitrogens, three methyls and two double bonds
:width: 420px
```

The same steps, if you'd rather read than watch:

| Step | Point at | Type | You get |
|---|---|---|---|
| 1 | empty canvas, with the **Rings** tool on cyclohexane | click | a cyclohexane |
| 2 | the ring's upper-right bond | `5` | a fused five-membered ring |
| 3 | the ring carbon at the lower right | `2` | a C=O |
| 4 | the carbon at the bottom left | `2` | the second C=O |
| 5 | the carbon between the two C=O | `n`, then `1` | an N with a methyl |
| 6 | the carbon at the left | `n`, then `1` | another N–CH₃ |
| 7 | the five-ring atom at the right | `n`, then `1` | the third N–CH₃ |
| 8 | the five-ring atom at the left | `n` | the last nitrogen |
| 9 | the bond shared by the two rings | `2` | a double bond |
| 10 | the bond from the top carbon to that nitrogen | `2` | the second double bond |

```{image} ../_static/tutorials/caffeine-keys.png
:alt: The finished caffeine
:width: 420px
```

Check it: the status bar reads C₈H₁₀N₄O₂, MW 194.19. A key on the wrong atom is one **Undo** (⌘Z or
Ctrl+Z) away.

The keys you used:

- On a bond: `5` fuses a ring of that size (`3`–`8` and `v` work too), `2` makes it double.
- On an atom: `2` adds a C=O to a ring carbon, `n` makes it nitrogen, `1` sprouts a single bond.

Every key is in [Keyboard shortcuts](../keys.md), and **Help → Keyboard Shortcuts** (F1) in the app.
Hover over a tool to see its own key.
