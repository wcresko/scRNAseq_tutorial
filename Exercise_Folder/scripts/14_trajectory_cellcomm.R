#!/usr/bin/env Rscript
# Talapas analysis pipeline 12 — parallels laptop Tutorial 12 (Trajectory & Cell-Cell Comm).
# Learning notebook: Exercise_Folder/Tutorial_14_Trajectory_CellCommunication.qmd
# Runs on the annotated ifnb object (the same dataset as the laptop notebook).
# Run:  sbatch --job-name=traj --mem=64G run_rscript.sbatch 14_trajectory_cellcomm.R
# In:   ../objects/ifnb_annotated.rds   Out: ../objects/ifnb_slingshot_pseudotime.csv
# Figures/tables: ../output/Mod14/Mod14_*  (Tutorial_14 has no executable chunks, so
# there are no canonical ModN_C{k}_ filenames to mirror — Mod14_<name> is used here.)

suppressPackageStartupMessages({
  library(Seurat); library(slingshot); library(SingleCellExperiment); library(tidyverse)
  library(patchwork)
})
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod14") # figures/tables for this module
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
seu <- readRDS(file.path(OBJ_DIR, "ifnb_annotated.rds"))
DefaultAssay(seu) <- "RNA"

# Part A — Pseudotime with Slingshot. The root flips the whole ordering, so you MUST
# justify it (progenitor markers / CytoTRACE / velocity); the placeholder below is
# only to demonstrate the mechanics.
sce  <- as.SingleCellExperiment(seu)
root <- levels(factor(seu$celltype_manual))[1]   # <- replace with a justified progenitor
sce  <- slingshot(sce, clusterLabels = "celltype_manual",
                  reducedDim = "UMAP", start.clus = root)
pt_tbl <- as.data.frame(slingPseudotime(sce)) |> rownames_to_column("cell")
write_csv(pt_tbl, file.path(OBJ_DIR, "ifnb_slingshot_pseudotime.csv"))
# Table out: same pseudotime table mirrored into the module output dir
write_csv(pt_tbl, file.path(OUT_DIR, "Mod14_slingshot_pseudotime.csv"))

# Figure out: Slingshot lineage curves over the UMAP, spots coloured by cell type.
# slingshot's curves live on the SingleCellExperiment already computed above, so this
# is one cheap base-graphics call (no new analysis). PNG via the base device.
umap <- reducedDims(sce)$UMAP
ct   <- factor(seu$celltype_manual)
pal  <- setNames(scales::hue_pal()(nlevels(ct)), levels(ct))
png(file.path(OUT_DIR, "Mod14_slingshot_umap_lineages.png"), width = 8, height = 6,
    units = "in", res = 300)
plot(umap, pch = 16, cex = 0.5, col = pal[ct],
     main = "Slingshot pseudotime lineages — ifnb (mechanics demo)",
     xlab = "UMAP 1", ylab = "UMAP 2")
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
legend("topright", legend = levels(ct), col = pal, pch = 16, cex = 0.6,
       title = "Cell type", bty = "n")
dev.off()

cat("Wrote Slingshot pseudotime for", ncol(sce), "cells (root =", root, ")\n")
cat("Wrote Mod14 figures/tables to", OUT_DIR, "\n")

# Part B — RNA velocity (scVelo) needs spliced/unspliced layers from velocyto/kb, NOT the
# filtered counts matrix. Python sketch — see the notebook:
#   scv.pp.filter_and_normalize(adata); scv.tl.velocity(adata, mode="dynamical"); ...

# Part C — Cell-cell communication with CellChat (uncomment once CellChat is installed):
# library(CellChat)
# cc <- createCellChat(seu, group.by = "celltype_manual"); cc@DB <- CellChatDB.human
# cc <- subsetData(cc); cc <- identifyOverExpressedGenes(cc)
# cc <- identifyOverExpressedInteractions(cc); cc <- computeCommunProb(cc)
# cc <- filterCommunication(cc, min.cells = 10); cc <- aggregateNet(cc)
# saveRDS(cc, file.path(OBJ_DIR, "ifnb_cellchat.rds"))
