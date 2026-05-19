#!/bin/bash
# Render script for the Single Cell RNA-seq Workshop site.
#
# Produces (by default) three kinds of output in docs/:
#   1. Lecture slides (revealjs)              — Lecture_Folder/_quarto.yml
#   2. Course website pages (html)            — root _quarto.yml
#      This includes the long-form reading material in Resources_Folder/
#      (appendices, Glossary, VS Code & Talapas) which renders as plain
#      HTML pages, not slides.
#   3. Tutorial .qmd sources copied verbatim into docs/Exercise_Folder/
#      so students can download them and run them locally. The student-
#      facing HTML for each tutorial is built by Exercise_Folder/_quarto.yml
#      with `execute.eval: false`, so the page shows the code but does NOT
#      run it. The Talapas Advanced sub-folder is rendered separately by
#      Exercise_Folder/Talapas/_quarto.yml.
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

echo "=== Rendering laptop tutorials (student version, eval=false) ==="
cd Exercise_Folder
quarto render
cd ..

echo "=== Rendering Talapas Advanced tutorials (eval=false) ==="
cd Exercise_Folder/Talapas
quarto render
cd ../..

# Copy each laptop-tutorial .qmd source to docs/Exercise_Folder/ as a *.qmd.txt
# file so students can download it cleanly. The .txt extension stops
# Quarto from auto-rewriting the link href to the rendered .html and
# stops the browser from trying to render the file inline; students are
# instructed (in Tutorials.qmd and each tutorial's intro) to drop the
# trailing `.txt` after downloading.
echo "=== Copying tutorial .qmd sources to docs/ as .qmd.txt for download ==="
mkdir -p docs/Exercise_Folder
for f in Exercise_Folder/*.qmd; do
  cp -f "$f" "docs/Exercise_Folder/$(basename "$f").txt"
done

# Talapas Advanced tutorial sources + executable scripts.
mkdir -p docs/Exercise_Folder/Talapas/scripts
for f in Exercise_Folder/Talapas/*.qmd; do
  [ -e "$f" ] || continue
  cp -f "$f" "docs/Exercise_Folder/Talapas/$(basename "$f").txt"
done
# Executable SLURM scripts and R scripts — copy as-is (not .txt-suffixed) so
# students get the right extensions when they download.
for f in Exercise_Folder/Talapas/scripts/*; do
  [ -e "$f" ] || continue
  cp -f "$f" "docs/Exercise_Folder/Talapas/scripts/$(basename "$f")"
done

# Companion: scATAC raw-data Cell Ranger tutorial (kept as-is in Talapas_RawData/).
mkdir -p docs/Exercise_Folder/Talapas_RawData
for f in Exercise_Folder/Talapas_RawData/*.qmd; do
  [ -e "$f" ] || continue
  cp -f "$f" "docs/Exercise_Folder/Talapas_RawData/$(basename "$f").txt"
done

# Copy metadata templates into the rendered docs so students can download them
mkdir -p docs/Resources_Folder/metadata_templates
cp -f Resources_Folder/metadata_templates/* docs/Resources_Folder/metadata_templates/ 2>/dev/null || true

if [ "${SOLUTIONS:-0}" = "1" ]; then
  echo "=== Rendering INSTRUCTOR SOLUTIONS (eval=true) ==="
  cd Exercise_Folder
  quarto render --profile solutions
  cd ..
  echo "Solutions written to docs/Exercise_Folder/_solutions/"
fi

echo "=== Done! ==="
echo "Preview: open docs/index.html in a browser."
