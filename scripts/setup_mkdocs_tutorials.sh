#!/usr/bin/env bash
# scripts/setup_mkdocs_tutorials.sh
# ─────────────────────────────────────────────────────────────────────────────
# Creates symlinks in docs_mkdocs/tutorials/ pointing to the notebooks stored
# under docs/src/notebooks/  (the Sphinx tutorial source of truth).
# Run once before `mkdocs serve` or `mkdocs build`.
# In CI this is replaced by a plain `cp` step in mkdocs.yaml.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$REPO_ROOT/docs/src/notebooks"
DST="$REPO_ROOT/docs_mkdocs/tutorials"

echo "Setting up MkDocs tutorial symlinks ..."
echo "  Source : $SRC"
echo "  Target : $DST"

mkdir -p "$DST"

for nb in "$SRC"/*.ipynb; do
    name="$(basename "$nb")"
    link="$DST/$name"
    if [ -L "$link" ]; then
        echo "  [skip]  $name (symlink already exists)"
    elif [ -e "$link" ]; then
        echo "  [skip]  $name (regular file already exists — remove it first)"
    else
        ln -s "$nb" "$link"
        echo "  [link]  $name"
    fi
done

# Symlink the img/ directory if present
if [ -d "$SRC/img" ] && [ ! -e "$DST/img" ]; then
    ln -s "$SRC/img" "$DST/img"
    echo "  [link]  img/"
fi

# Symlink the logo into docs_mkdocs/assets/
ASSETS="$REPO_ROOT/docs_mkdocs/assets"
mkdir -p "$ASSETS"
LOGO="$REPO_ROOT/docs/src/neko_logo.png"
if [ -f "$LOGO" ] && [ ! -e "$ASSETS/neko_logo.png" ]; then
    ln -s "$LOGO" "$ASSETS/neko_logo.png"
    echo "  [link]  assets/neko_logo.png"
fi

echo "Done. You can now run:  mkdocs serve"
