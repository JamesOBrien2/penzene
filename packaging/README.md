# Packaging

Files prepared for distribution channels outside this repository. Each is submitted to its
channel by hand, with the owner's approval.

| Channel | File | Submitted? |
|---|---|---|
| conda-forge | `conda-forge/recipe.yaml`: the Python module, built against conda-forge's RDKit and Qt | not yet |
| Homebrew | `homebrew/penzene.rb.in`: a cask for the macOS disk images ([README](homebrew/README.md)) | not yet; needs notarization (#49) |
| winget | `winget/*.yaml.in`: manifests for the Windows installer ([README](winget/README.md)) | not yet; ideally after signing (#110) |

Test a recipe locally with rattler-build (a variants file gives what conda-forge's pinning
normally supplies, e.g. `c_stdlib` on macOS):

```sh
pixi exec rattler-build build -r packaging/conda-forge/recipe.yaml -c conda-forge -m variants.yaml
```

Code signing for the macOS and Windows release packages, the secrets it needs and the
Windows alternatives (Azure Trusted Signing, SignPath): [SIGNING.md](SIGNING.md).
