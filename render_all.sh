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
#      run it. All tutorials — core 00-10 (incl. the Talapas HPC modules)
#      and bonus 11-17 — now live flat in Exercise_Folder/, with the
#      executable Talapas scripts in Exercise_Folder/scripts/.
#
# Optionally produce instructor "solutions" by setting SOLUTIONS=1 in the
# environment. That re-renders the tutorials with `eval: true` (so chunks
# actually run) into docs/Exercise_Folder/_solutions/. Requires the
# workshop dataset under ../data/ and all packages from Software_Setup.qmd.
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

echo "=== Rendering tutorials — core 00-10 + bonus 11-17 (student version, eval=false) ==="
cd Exercise_Folder
quarto render
cd ..

# Copy each tutorial .qmd source to docs/Exercise_Folder/ as a *.qmd.txt
# file so students can download it cleanly. The .txt extension stops
# Quarto from auto-rewriting the link href to the rendered .html and
# stops the browser from trying to render the file inline; students are
# instructed (in Materials.qmd and each tutorial's intro) to drop the
# trailing `.txt` after downloading. This covers every tutorial — the
# core analysis modules (00-08), the Talapas HPC modules (09-10), and the
# bonus modules (11-17) — since they all now live flat in Exercise_Folder/.
echo "=== Copying tutorial .qmd sources to docs/ as .qmd.txt for download ==="
mkdir -p docs/Exercise_Folder
for f in Exercise_Folder/*.qmd; do
  cp -f "$f" "docs/Exercise_Folder/$(basename "$f").txt"
done

# Executable SLURM (.sbatch) and R (.R) scripts — the Talapas track's fourth
# leg. Copied as-is (not .txt-suffixed) so students get the right extensions.
echo "=== Copying Talapas scripts to docs/ ==="
mkdir -p docs/Exercise_Folder/scripts
for f in Exercise_Folder/scripts/*; do
  [ -e "$f" ] || continue
  cp -f "$f" "docs/Exercise_Folder/scripts/$(basename "$f")"
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
