#!/usr/bin/env Rscript
# Talapas analysis pipeline 14 — parallels laptop Tutorial 14 (Trajectory & Cell-Cell Comm).
# RNA velocity (scVelo) is Python — see the companion script 14_rna_velocity_scvelo.py.
# Learning notebook: Exercise_Folder/Tutorial_14_Trajectory_CellCommunication.qmd
# Runs on the annotated ifnb object (the same dataset as the laptop notebook).
# Run:  sbatch --job-name=traj --mem=64G run_rscript.sbatch 14_trajectory_cellcomm.R
# In:   ../data/ifnb_annotated.rds   Out: ../data/ifnb_slingshot_pseudotime.csv
# Figures/tables: ../output/Mod14/Mod14_*  (Tutorial_14 has no executable chunks, so
# there are no canonical ModN_C{k}_ filenames to mirror — Mod14_<name> is used here.)

suppressPackageStartupMessages({
  library(Seurat); library(slingshot); library(SingleCellExperiment); library(tidyverse)
  library(patchwork)
})
set.seed(2026)
DATA_DIR <- Sys.getenv("DATA_DIR", "../data")
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod14") # figures/tables for this module
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)   # pipeline hand-off objects (.rds/.csv)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("[dirs] data -> ", normalizePath(DATA_DIR), "  |  figures/tables -> ", normalizePath(OUT_DIR))

# save_base_fig() does the same for BASE-graphics figures (plot()/heatmap()/etc.):
# pass a no-argument function that draws the figure; it is rendered to .png and .svg.
save_base_fig <- function(filename, draw, width = 7, height = 5, res = 300) {
  grDevices::png(filename, width = width, height = height, units = "in", res = res)
  draw(); grDevices::dev.off()
  svg_path <- paste0(tools::file_path_sans_ext(filename), ".svg")
  tryCatch({ grDevices::svg(svg_path, width = width, height = height); draw(); grDevices::dev.off() },
           error = function(e) message("  (could not write ", basename(svg_path), ")"))
}
seu <- readRDS(file.path(DATA_DIR, "ifnb_annotated.rds"))
DefaultAssay(seu) <- "RNA"

# Part A — Pseudotime with Slingshot. The root flips the whole ordering, so you MUST
# justify it (progenitor markers / CytoTRACE / velocity); the placeholder below is
# only to demonstrate the mechanics.
sce  <- as.SingleCellExperiment(seu)
root <- levels(factor(seu$celltype_manual))[1]   # <- replace with a justified progenitor
sce  <- slingshot(sce, clusterLabels = "celltype_manual",
                  reducedDim = "UMAP", start.clus = root)
pt_tbl <- as.data.frame(slingPseudotime(sce)) |> rownames_to_column("cell")
write_csv(pt_tbl, file.path(DATA_DIR, "ifnb_slingshot_pseudotime.csv"))
# Table out: same pseudotime table mirrored into the module output dir
write_csv(pt_tbl, file.path(OUT_DIR, "Mod14_slingshot_pseudotime.csv"))

# Figure out: Slingshot lineage curves over the UMAP, spots coloured by cell type.
# slingshot's curves live on the SingleCellExperiment already computed above, so this
# is one cheap base-graphics call (no new analysis). PNG via the base device.
umap <- reducedDims(sce)$UMAP
ct   <- factor(seu$celltype_manual)
pal  <- setNames(scales::hue_pal()(nlevels(ct)), levels(ct))
save_base_fig(file.path(OUT_DIR, "Mod14_slingshot_umap_lineages.png"), width = 8, height = 6, draw = function() {
  plot(umap, pch = 16, cex = 0.5, col = pal[ct],
       main = "Slingshot pseudotime lineages — ifnb (mechanics demo)",
       xlab = "UMAP 1", ylab = "UMAP 2")
  lines(SlingshotDataSet(sce), lwd = 2, col = "black")
  legend("topright", legend = levels(ct), col = pal, pch = 16, cex = 0.6,
         title = "Cell type", bty = "n")
})

cat("Wrote Slingshot pseudotime for", ncol(sce), "cells (root =", root, ")\n")
cat("Wrote Mod14 figures/tables to", OUT_DIR, "\n")

# Part B — RNA velocity (scVelo) is a Python analysis. Run the companion script
# 14_rna_velocity_scvelo.py (on Talapas: sbatch run_python.sbatch
# 14_rna_velocity_scvelo.py). It needs spliced/unspliced layers from velocyto/kb,
# NOT the Cell Ranger filtered matrix.

# Part C — Cell-cell communication with CellChat (uncomment once CellChat is installed):
# library(CellChat)
# cc <- createCellChat(seu, group.by = "celltype_manual"); cc@DB <- CellChatDB.human
# cc <- subsetData(cc); cc <- identifyOverExpressedGenes(cc)
# cc <- identifyOverExpressedInteractions(cc); cc <- computeCommunProb(cc)
# cc <- filterCommunication(cc, min.cells = 10); cc <- aggregateNet(cc)
# saveRDS(cc, file.path(DATA_DIR, "ifnb_cellchat.rds"))
