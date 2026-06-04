#!/usr/bin/env Rscript
# Talapas analysis pipeline 02 — parallels laptop Tutorial 02 (DimReduction & Clustering).
# Learning notebook: Exercise_Folder/Tutorial_02_DimReduction_Clustering.qmd
# Run:  sbatch --job-name=cluster run_rscript.sbatch 02_dimreduction_clustering.R
# In:   ../objects/ifnb_preprocessed.rds   Out: ../objects/ifnb_clustered.rds

suppressPackageStartupMessages({ library(Seurat); library(tidyverse) })
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "ifnb_preprocessed.rds"))

# Step 1 — PCA
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

# Step 2 — Choose PCs (inspect ElbowPlot interactively; 30 is a safe default here)
N_PCS <- 30

# Step 3 — Neighbor graph
seu <- FindNeighbors(seu, dims = 1:N_PCS, verbose = FALSE)

# Step 4 — Cluster at several resolutions
# FindClusters() uses the Louvain default (algorithm = 1); add algorithm = 4 for
# Leiden (needs the leidenalg Python package via reticulate). See Lecture 02.
for (res in c(0.1, 0.3, 0.5, 0.7, 1.0))
  seu <- FindClusters(seu, resolution = res, verbose = FALSE)
Idents(seu) <- "RNA_snn_res.0.5"

# Step 5 — UMAP
seu <- RunUMAP(seu, dims = 1:N_PCS, verbose = FALSE)

saveRDS(seu, file.path(OBJ_DIR, "ifnb_clustered.rds"))
cat("Wrote", file.path(OBJ_DIR, "ifnb_clustered.rds"), "with",
    length(levels(Idents(seu))), "clusters\n")
