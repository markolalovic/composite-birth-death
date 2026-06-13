#!/usr/bin/env bash
set -euo pipefail

# render Quarto docs
quarto render docs-src

# disable jekyll processing on github pages
touch docs/.nojekyll

echo "Site built in docs/"