#!/bin/bash
# Publish the repository on GitHub and create the release that Zenodo archives.
#
# PREREQUISITE: the empty repository must already exist at
#   https://github.com/wellfonseca/schrodinger-fem
# Create it at https://github.com/new with no README, no .gitignore and no
# license, so that it does not conflict with what is already committed here.
#
# Usage:
#   ./PUBLISH.sh 1.0.0
set -euo pipefail
cd "$(dirname "$0")"
VER="${1:-1.0.0}"
REMOTE="git@github.com:wellfonseca/schrodinger-fem.git"

echo "==> remote"
git remote remove origin 2>/dev/null || true
git remote add origin "$REMOTE"

echo "==> pushing main"
git push -u origin main

echo "==> creating the annotated tag v$VER"
git tag -a "v$VER" -m "schrodinger-fem v$VER"
git push origin "v$VER"

echo
echo "Tag v$VER pushed."
echo
echo "NEXT STEP (manual, in the browser):"
echo "  1. https://github.com/wellfonseca/schrodinger-fem/releases/new"
echo "     - Tag: v$VER"
echo "     - Title: schrodinger-fem v$VER"
echo "     - Description: paste RELEASE_NOTES_v$VER.md"
echo "     - Publish release"
echo "  2. Zenodo mints the DOI automatically within a few minutes."
echo "  3. Check https://zenodo.org/account/settings/github/"
echo "  4. Put the DOI in the manuscript, section Code availability."
