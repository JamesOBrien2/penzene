#!/usr/bin/env bash
# Installs a self-contained Penzene.app in /Applications, for trying a build out locally:
#
#   pixi run install-app
#
# Replaces an existing Penzene.app there, and removes the ~/Applications copy that earlier
# versions of this task installed. PENZENE_INSTALL_DIR installs somewhere else instead.
set -euo pipefail
dest=${PENZENE_INSTALL_DIR:-/Applications}
app="$dest/Penzene.app"
staged="$dest/.Penzene.app.installing"
replaced="$dest/.Penzene.app.replaced"
old="$HOME/Applications/Penzene.app"

fail() {
    rm -rf "$staged" 2>/dev/null || true
    echo "error: can't install $app: $1" >&2
    echo "Installing in $dest needs write access there: use an administrator account, or" >&2
    echo "set PENZENE_INSTALL_DIR to a folder you can write to." >&2
    exit 1
}

rm -rf build/deploy && mkdir -p build/deploy && cp -R build/bin/penzene.app build/deploy/
macdeployqt6 build/deploy/penzene.app -libpath="$CONDA_PREFIX/lib" 2>/dev/null || true
codesign --force --deep --sign - build/deploy/penzene.app

# Copy next to the old app and swap by renaming, so a failure leaves the old one working.
[ -d "$dest" ] && [ -w "$dest" ] || fail "$dest isn't writable"
rm -rf "$staged" "$replaced" 2>/dev/null || true
cp -R build/deploy/penzene.app "$staged" || fail "copying failed"
if [ -e "$app" ]; then mv "$app" "$replaced" || fail "the existing app couldn't be replaced"; fi
mv "$staged" "$app" || { mv "$replaced" "$app" 2>/dev/null; fail "moving the new app into place failed"; }
echo "Installed $app"
rm -rf "$replaced" 2>/dev/null || echo "warning: couldn't delete the old app, left at $replaced" >&2

if [ -e "$old" ] && [ "$(cd "$dest" && pwd -P)" != "$(cd "$HOME/Applications" && pwd -P)" ]; then
    rm -rf "$old" && echo "Removed the old copy in ~/Applications"
fi
