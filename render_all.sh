#!/bin/bash
# Render script for the Single Cell RNA-seq Workshop site.
#
# Produces three kinds of output in docs/:
#   1. Lecture slides (revealjs)      — Lecture_Folder/_quarto.yml
#   2. Course website pages (html)    — root _quarto.yml
#   3. Exercise HTML + downloadable .qmd sources side-by-side:
#        docs/Exercise_Folder/Tutorial_01_QC_Preprocessing.html
#        docs/Exercise_Folder/Tutorial_01_QC_Preprocessing.qmd  (download)
#
# Exercises are rendered with `execute.eval: false` in their YAML, so the
# HTML shows the code but does NOT run it. The .qmd files are copied
# verbatim via the `resources:` entry in the root _quarto.yml so students
# can download and open them locally in RStudio to actually run the code.
#
# Usage: bash render_all.sh

set -e

echo "=== Rendering lecture slides (revealjs) ==="
cd Lecture_Folder
quarto render
cd ..

echo "=== Rendering course website (pages + exercises) ==="
quarto render

# Belt-and-braces: in case the `resources:` copy in _quarto.yml does not
# propagate on some Quarto versions, make sure every exercise .qmd ends
# up next to its rendered .html for students to download.
echo "=== Ensuring exercise .qmd sources are available for download ==="
mkdir -p docs/Exercise_Folder
cp -f Exercise_Folder/*.qmd docs/Exercise_Folder/

echo "=== Done! ==="
echo "Preview: open docs/index.html in a browser."
