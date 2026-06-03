#!/usr/bin/env Rscript
# Talapas 07 — Multi-Sample Integration with Harmony
# Treats each donor (from CMO multiplexing) as a separate "sample" and
# corrects donor-level batch effects on the PCA embedding using harmony.
# Writes nsclc_integrated.rds.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(tidyverse)
})

set.seed(2026)

OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "nsclc_annotated.rds"))

# If the object lacks per-cell donor IDs (e.g. CMO assignment was not run),
# synthesize 7 pseudo-donors to demonstrate the workflow.
if (!"donor" %in% colnames(seu@meta.data)) {
  seu$donor <- paste0("d", sample(1:7, ncol(seu), replace = TRUE))
}

# Re-run PCA on the (already-normalized) joint matrix
DefaultAssay(seu) <- "RNA"
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

# Harmony correction over the donor variable
seu <- RunHarmony(seu, group.by.vars = "donor",
                  reduction.use = "pca", reduction.save = "harmony",
                  verbose = FALSE)

# Re-cluster + UMAP on the harmony embedding
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30,
               reduction.name = "umap_harmony", verbose = FALSE)

saveRDS(seu, file.path(OBJ_DIR, "nsclc_integrated.rds"))
cat("Wrote", file.path(OBJ_DIR, "nsclc_integrated.rds"), "\n")
