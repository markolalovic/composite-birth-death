#!/usr/bin/env bash
set -euo pipefail

# render Quarto docs
quarto render docs-src

# add google search verification file
cp docs-src/google204cb8c5f11feeb9.html docs/

# disable jekyll processing on github pages
touch docs/.nojekyll

echo "Site built in docs/"
