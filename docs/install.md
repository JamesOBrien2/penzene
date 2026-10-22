# Installing

Choose your system to download the latest stable installer directly:

| System | File |
|---|---|
| macOS, Apple silicon | [Download DMG](https://github.com/JamesOBrien2/penzene/releases/latest/download/penzene-macos-arm64.dmg) |
| macOS, Intel | [Download DMG](https://github.com/JamesOBrien2/penzene/releases/latest/download/penzene-macos-x86_64.dmg) |
| Windows | [Download installer](https://github.com/JamesOBrien2/penzene/releases/latest/download/penzene-windows-x64-setup.exe) or [portable ZIP](https://github.com/JamesOBrien2/penzene/releases/latest/download/penzene-windows-x64.zip) |
| Linux | [Download AppImage](https://github.com/JamesOBrien2/penzene/releases/latest/download/penzene-linux-x86_64.AppImage) |

These links always point to the newest stable release. The installers are free to use and
Penzene is licensed under [GPL-3.0-or-later](https://github.com/JamesOBrien2/penzene/blob/main/LICENSE).
The distributed dependencies and their notices are listed in
[`THIRD_PARTY_LICENSES.md`](https://github.com/JamesOBrien2/penzene/blob/main/THIRD_PARTY_LICENSES.md).
Universities may deploy the installers through their usual software distribution tools; consult
your institution's software licensing process when redistributing copies.

On a Mac, **Apple menu → About This Mac** tells you whether its chip is Apple silicon (M-series)
or Intel. On Windows, choose the **installer** for ordinary use; the portable ZIP is for a
computer where you cannot install software and does not register Penzene's editable drawing
object in Word or PowerPoint.

## Opening it the first time

:::{note}
Penzene's installers aren't signed with a paid Apple or Microsoft certificate, so each system
asks once whether to trust it.
:::

**macOS.** Drag Penzene to Applications and open it. If macOS says it "can't be opened" or can't be
checked, go to **System Settings → Privacy & Security**, scroll down to the message about Penzene,
and click **Open Anyway** after confirming you downloaded it from the release linked above. You
only need to do this once. On a managed university computer, ask your IT team if that control is
unavailable.

**Windows.** If SmartScreen shows "Windows protected your PC", confirm the installer came from the
release linked above, then click **More info → Run anyway**.

**Linux.** Make the AppImage executable, then run it:

```sh
chmod +x penzene-linux-x86_64.AppImage
./penzene-linux-x86_64.AppImage
```

## Python package

```sh
pip install penzene
```

Python 3.12+ on macOS, Linux (glibc 2.17+) and Windows. See [Python](python.md).

## Building from source

See [Contributing](contributing.md).
