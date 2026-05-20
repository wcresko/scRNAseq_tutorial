#!/usr/bin/env Rscript
# Talapas 04 — Dim Reduction & Clustering
# Reads nsclc_preprocessed.rds, runs PCA, neighbor graph, Leiden clustering,
# and UMAP. Writes nsclc_clustered.rds for the next step.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

set.seed(2026)

OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "nsclc_preprocessed.rds"))

# PCA
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

# Pick PCs (default 30; adjust after inspecting ElbowPlot interactively)
N_PCS <- 30

# Neighbor graph + clustering at multiple resolutions
seu <- FindNeighbors(seu, dims = 1:N_PCS, verbose = FALSE)
for (res in c(0.3, 0.5, 0.8)) {
  seu <- FindClusters(seu, resolution = res, verbose = FALSE)
}
Idents(seu) <- "RNA_snn_res.0.5"

# UMAP
seu <- RunUMAP(seu, dims = 1:N_PCS, verbose = FALSE)

saveRDS(seu, file.path(OBJ_DIR, "nsclc_clustered.rds"))
cat("Wrote", file.path(OBJ_DIR, "nsclc_clustered.rds"),
    "with", length(levels(Idents(seu))), "clusters\n")
