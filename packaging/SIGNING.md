# Code signing

`.github/workflows/release.yml` signs the macOS and Windows packages when the secrets below
exist (Settings → Secrets and variables → Actions → New repository secret). Without them the
signing steps are skipped and the release builds exactly as before: the macOS app ad-hoc
signed, the Windows files unsigned. Pull requests from forks never see secrets, so they build
unsigned too.

## macOS: Developer ID and notarization (#49)

**Cost:** the Apple Developer Program, US$99 a year (developer.apple.com/programs). The
Account Holder must create the Developer ID certificate.

| Secret | What it is |
|---|---|
| `MACOS_CERTIFICATE` | the Developer ID Application certificate and its private key, as a base64 .p12 |
| `MACOS_CERTIFICATE_PASSWORD` | the password chosen when exporting the .p12 |
| `NOTARY_KEY` | an App Store Connect API key, the base64 of its `AuthKey_XXXXXXXXXX.p8` |
| `NOTARY_KEY_ID` | that key's Key ID (10 characters) |
| `NOTARY_ISSUER_ID` | the Issuer ID (a UUID) shown above the key list |

Signing needs only the first two; notarization also needs the three `NOTARY_*` secrets.

1. **Certificate.** In Xcode → Settings → Accounts → Manage Certificates, add a *Developer ID
   Application* certificate (or create one at developer.apple.com → Certificates with a CSR
   from Keychain Access). In Keychain Access → My Certificates, right-click
   "Developer ID Application: …" → Export as .p12 with a password, then:
   ```sh
   base64 -i DeveloperID.p12 | pbcopy   # paste as MACOS_CERTIFICATE
   ```
2. **API key.** App Store Connect → Users and Access → Integrations → App Store Connect API →
   Team Keys → generate a key with the *Developer* role. Download the .p8 (possible once only),
   note the Key ID and the Issuer ID, then:
   ```sh
   base64 -i AuthKey_XXXXXXXXXX.p8 | pbcopy   # paste as NOTARY_KEY
   ```

The workflow imports the certificate into a temporary keychain, signs `penzene.app` with the
hardened runtime (`codesign --deep --force --options runtime --timestamp`, no entitlements,
since nothing in the app needs JIT or foreign libraries), smoke-tests the signed app, builds
the dmg, signs it, submits it with `xcrun notarytool submit --wait` (printing Apple's log if
it is rejected) and staples the ticket.

**Verify** a downloaded dmg:

```sh
codesign --verify --deep --strict --verbose=2 /Volumes/Penzene/penzene.app
spctl -a -vv /Volumes/Penzene/penzene.app          # "source=Notarized Developer ID"
spctl -a -vv -t open --context context:primary-signature penzene-macos-arm64.dmg
xcrun stapler validate penzene-macos-arm64.dmg
```

## Windows: Authenticode (#110)

The workflow signs `penzene.exe` (before it is zipped and packed into the installer) and
`penzene-windows-x64-setup.exe` with `signtool`, SHA-256 and the DigiCert RFC 3161 timestamp
server.

| Secret | What it is |
|---|---|
| `WINDOWS_CERTIFICATE` | a code-signing certificate and private key, as a base64 .pfx |
| `WINDOWS_CERTIFICATE_PASSWORD` | the .pfx password |

Export and encode it (PowerShell):

```powershell
[Convert]::ToBase64String([IO.File]::ReadAllBytes("penzene.pfx")) | Set-Clipboard
```

**Costs and the catch.** Since June 2023 certificate authorities issue OV and EV code-signing
keys only on hardware (a USB token or the CA's cloud HSM), so a new certificate usually
cannot be exported as a .pfx. The .pfx path suits a certificate that can be exported (an
older one, or a CA offering one). Otherwise:

- **Azure Trusted Signing** (recommended): about US$10 a month (Basic tier), Microsoft-managed
  keys, reputation builds up for SmartScreen. Needs an Azure subscription and identity
  validation (an organisation, or an individual in the US or Canada at present). It would
  replace the two "Sign" steps with `azure/trusted-signing-action`, authenticated with an
  Azure app registration (secrets `AZURE_TENANT_ID`, `AZURE_CLIENT_ID`,
  `AZURE_CLIENT_SECRET`, plus the account's endpoint, account name and certificate profile).
  Not implemented here; say so if this is the route taken.
- **SignPath.io**: free for open-source projects that qualify, signing through their service.
- **OV/EV certificate from a CA** on a cloud HSM (e.g. DigiCert KeyLocker, SSL.com eSigner):
  roughly US$200–600 a year, each with its own signing tool in place of `signtool /f`.

An OV certificate still shows SmartScreen warnings until the file earns reputation; EV and
Trusted Signing build it faster.

**Verify** on Windows:

```powershell
signtool verify /pa /v penzene-windows-x64-setup.exe
signtool verify /pa /v "C:\Program Files\Penzene\penzene.exe"
```

or right-click the file → Properties → Digital Signatures.

The Inno Setup uninstaller (`unins000.exe`) stays unsigned; sign it too by moving the
installer signing into `penzene.iss` (`SignTool=` with `iscc /S`) if that matters.
