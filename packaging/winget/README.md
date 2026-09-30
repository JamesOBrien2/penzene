# winget

The three `*.yaml.in` files are the winget manifests (schema 1.12.0) for the Inno Setup installer,
`penzene-windows-x64-setup.exe`. `render.py` fills in the version, release date and the
installer's sha256 (from GitHub's asset digest, or by downloading it), in the directory layout
microsoft/winget-pkgs uses:

```sh
python3 packaging/winget/render.py 1.4.0 out   # writes out/manifests/j/JamesOBrien2/Penzene/1.4.0/
```

Installer details, from `cmake/penzene.iss`: there's no explicit `AppId`, so Inno registers the
uninstaller as `Penzene_is1` (the `ProductCode`). The installer can run per user
(`/CURRENTUSER`) or for all users (`/ALLUSERS`), so there's one entry for each scope; the silent
switches are `/VERYSILENT /SUPPRESSMSGBOXES`.

**Ideally after signing (#110).** winget accepts unsigned installers, but SmartScreen and
Defender warn about them more often.

## Submitting

1. Test them: the `winget` workflow (`.github/workflows/winget.yml`, run by hand from the Actions
   tab for a release) renders them, runs `winget validate`, installs from them and uninstalls on
   Windows, the two checks the winget-pkgs pull request checklist asks for. On Windows the same by
   hand (after `winget settings --enable LocalManifestFiles` once, as administrator):

   ```powershell
   winget validate --manifest out\manifests\j\JamesOBrien2\Penzene\1.4.0
   winget install --manifest out\manifests\j\JamesOBrien2\Penzene\1.4.0
   winget uninstall --manifest out\manifests\j\JamesOBrien2\Penzene\1.4.0
   ```

2. Submit, either:
   - `wingetcreate submit out\manifests\j\JamesOBrien2\Penzene\1.4.0` (it asks for a GitHub
     token, forks microsoft/winget-pkgs and opens the pull request), or
   - by hand: fork [microsoft/winget-pkgs](https://github.com/microsoft/winget-pkgs), copy the
     `manifests/` tree into it, and open a pull request.

   Title it `New package: JamesOBrien2.Penzene version X.Y.Z` and tick the template's checklist
   (link the `winget` workflow run for validate and install). The bot closes a pull request that
   waits on its author: sign the Microsoft CLA when it asks, and answer any `Needs-Author-Feedback`
   within five days.
3. For later releases, `wingetcreate update JamesOBrien2.Penzene --version X.Y.Z --urls
   <installer url> --submit` updates the existing manifests, or render again and repeat step 2.
   Then add `winget install JamesOBrien2.Penzene` to `docs/install.md`.
