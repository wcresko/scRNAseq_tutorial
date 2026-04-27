#!/bin/bash
# Render script for the Single Cell RNA-seq Workshop site.
#
# Produces (by default) three kinds of output in docs/:
#   1. Lecture slides (revealjs)              — Lecture_Folder/_quarto.yml
#   2. Course website pages (html)            — root _quarto.yml
#      This includes the long-form reading material in Resources_Folder/
#      (appendices, Glossary, VS Code & Talapas) which renders as plain
#      HTML pages, not slides.
#   3. Exercise .qmd sources copied verbatim into docs/Exercise_Folder/
#      so students can download them and run them locally. The student-
#      facing HTML for each exercise is built by Exercise_Folder/_quarto.yml
#      with `execute.eval: false`, so the page shows the code but does NOT
#      run it.
#
# Optionally produce instructor "solutions" by setting SOLUTIONS=1 in the
# environment. That re-renders the tutorials with `eval: true` (so chunks
# actually run) into docs/Exercise_Folder/_solutions/. Requires the
# workshop dataset under ../data/ and all packages from Software.qmd.
#
# Usage:
#     bash render_all.sh                # student build
#     SOLUTIONS=1 bash render_all.sh    # student build + instructor solutions

set -e

echo "=== Rendering lecture slides (revealjs) ==="
cd Lecture_Folder
quarto render
cd ..

echo "=== Rendering course website (pages + Resources_Folder) ==="
quarto render

echo "=== Rendering exercise tutorials (student version, eval=false) ==="
cd Exercise_Folder
quarto render
cd ..

# Belt-and-braces: in case the `resources:` copy in root _quarto.yml does
# not propagate on some Quarto versions, make sure every exercise .qmd
# also ends up next to its rendered .html for students to download.
echo "=== Ensuring exercise .qmd sources are available for download ==="
mkdir -p docs/Exercise_Folder
cp -f Exercise_Folder/*.qmd docs/Exercise_Folder/

if [ "${SOLUTIONS:-0}" = "1" ]; then
  echo "=== Rendering INSTRUCTOR SOLUTIONS (eval=true) ==="
  cd Exercise_Folder
  quarto render --profile solutions
  cd ..
  echo "Solutions written to docs/Exercise_Folder/_solutions/"
fi

echo "=== Done! ==="
echo "Preview: open docs/index.html in a browser."
