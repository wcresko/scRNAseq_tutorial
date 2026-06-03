#!/usr/bin/env Rscript
# Talapas analysis pipeline 12 — parallels laptop Tutorial 12 (Trajectory & Cell-Cell Comm).
# Learning notebook: Exercise_Folder/Tutorial_12_Trajectory_CellCommunication.qmd
# Runs on the annotated NSCLC object (the laptop notebook uses the annotated ifnb).
# Run:  sbatch --job-name=traj --mem=64G run_rscript.sbatch 12_trajectory_cellcomm.R
# In:   ../objects/nsclc_annotated.rds   Out: ../objects/nsclc_slingshot_pseudotime.csv

suppressPackageStartupMessages({
  library(Seurat); library(slingshot); library(SingleCellExperiment); library(tidyverse)
})
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "nsclc_annotated.rds"))
DefaultAssay(seu) <- "RNA"

# Part A — Pseudotime with Slingshot. The root flips the whole ordering, so you MUST
# justify it (progenitor markers / CytoTRACE / velocity); the placeholder below is
# only to demonstrate the mechanics.
sce  <- as.SingleCellExperiment(seu)
root <- levels(factor(seu$celltype_manual))[1]   # <- replace with a justified progenitor
sce  <- slingshot(sce, clusterLabels = "celltype_manual",
                  reducedDim = "UMAP", start.clus = root)
write_csv(as.data.frame(slingPseudotime(sce)) |> rownames_to_column("cell"),
          file.path(OBJ_DIR, "nsclc_slingshot_pseudotime.csv"))
cat("Wrote Slingshot pseudotime for", ncol(sce), "cells (root =", root, ")\n")

# Part B — RNA velocity (scVelo) needs spliced/unspliced layers from velocyto/kb, NOT the
# filtered counts matrix. Python sketch — see the notebook:
#   scv.pp.filter_and_normalize(adata); scv.tl.velocity(adata, mode="dynamical"); ...

# Part C — Cell-cell communication with CellChat (uncomment once CellChat is installed):
# library(CellChat)
# cc <- createCellChat(seu, group.by = "celltype_manual"); cc@DB <- CellChatDB.human
# cc <- subsetData(cc); cc <- identifyOverExpressedGenes(cc)
# cc <- identifyOverExpressedInteractions(cc); cc <- computeCommunProb(cc)
# cc <- filterCommunication(cc, min.cells = 10); cc <- aggregateNet(cc)
# saveRDS(cc, file.path(OBJ_DIR, "nsclc_cellchat.rds"))
