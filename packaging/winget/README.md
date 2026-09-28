# winget

The three `*.yaml.in` files are the winget manifests (schema 1.6.0) for the Inno Setup installer,
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

1. On Windows, render the manifests and test them locally (needs `winget settings --enable
   LocalManifestFiles` once, as administrator):

   ```powershell
   winget validate --manifest out\manifests\j\JamesOBrien2\Penzene\1.4.0
   winget install --manifest out\manifests\j\JamesOBrien2\Penzene\1.4.0
   winget uninstall --id JamesOBrien2.Penzene
   ```

2. Submit, either:
   - `wingetcreate submit out\manifests\j\JamesOBrien2\Penzene\1.4.0` (it asks for a GitHub
     token, forks microsoft/winget-pkgs and opens the pull request), or
   - by hand: fork [microsoft/winget-pkgs](https://github.com/microsoft/winget-pkgs), copy the
     `manifests/` tree into it, and open a pull request.
3. For later releases, `wingetcreate update JamesOBrien2.Penzene --version X.Y.Z --urls
   <installer url> --submit` updates the existing manifests, or render again and repeat step 2.
   Then add `winget install JamesOBrien2.Penzene` to `docs/install.md`.
