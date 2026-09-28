# Homebrew cask

`penzene.rb.in` is the cask for the two macOS disk images; `render.py` fills in the version and
the images' sha256 (from GitHub's asset digests, or by downloading them).

```sh
python3 packaging/homebrew/render.py 1.4.0 > penzene.rb
```

**Blocked on signing (#49).** homebrew-cask no longer accepts apps that aren't signed with a
Developer ID and notarized, so submit once a notarized release is out.

## Submitting

1. Render the cask for the notarized release and check it installs, opens and uninstalls cleanly
   (`--no-quarantine` isn't needed once it is notarized):

   ```sh
   brew tap-new $USER/local && cp penzene.rb "$(brew --repo $USER/local)/Casks/penzene.rb"
   brew install --cask $USER/local/penzene && open -a penzene
   brew uninstall --cask --zap $USER/local/penzene
   brew audit --new --cask $USER/local/penzene && brew style --cask $USER/local/penzene
   ```

2. Fork [Homebrew/homebrew-cask](https://github.com/Homebrew/homebrew-cask), add the file as
   `Casks/p/penzene.rb`, and open a pull request following its template.
3. After it merges, Homebrew's autobump updates the version and checksums for new releases (the
   `livecheck` block follows the latest GitHub release). Add `brew install --cask penzene` to
   `docs/install.md`.

The app in the disk image is `penzene.app` (lower case); `zap` removes its settings
(`com.penzene.Penzene.plist`, from Qt's organisation and app names), its support directory
(`~/Library/Application Support/Penzene`) and saved window state.
