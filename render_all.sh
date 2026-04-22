#!/bin/bash
# Render script for BioE Stats website + lecture slides
#
# The lectures are a separate Quarto project (Lecture_Folder/_quarto.yml)
# so they render as true revealjs presentations, not website HTML pages.
#
# Usage: bash render_all.sh

set -e

echo "=== Rendering lecture slides (revealjs) ==="
cd Lecture_Folder
quarto render
cd ..

echo "=== Rendering course website ==="
quarto render

echo "=== Done! ==="
