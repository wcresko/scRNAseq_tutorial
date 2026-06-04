#!/usr/bin/env Rscript
# Talapas analysis pipeline 05 — parallels laptop Tutorial 05 (Integration, Harmony).
# Learning notebook: Exercise_Folder/Tutorial_05_Integration.qmd
# Run:  sbatch --job-name=integrate run_rscript.sbatch 05_integration.R
# In:   ../objects/ifnb_annotated_final.rds   Out: ../objects/ifnb_integrated.rds

suppressPackageStartupMessages({ library(Seurat); library(harmony); library(tidyverse) })
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
# Prefer the reference-reconciled object from step 04; fall back to step 03.
in_file <- file.path(OBJ_DIR, "ifnb_annotated_final.rds")
if (!file.exists(in_file)) in_file <- file.path(OBJ_DIR, "ifnb_annotated.rds")
seu <- readRDS(in_file)
DefaultAssay(seu) <- "RNA"

# ifnb's batch / sample axis is `stim` (CTRL vs IFN-beta STIM) — the two samples we
# integrate. `donor` (the 8 real `ind` patients) is carried for downstream DE / DA.
if (!"donor" %in% colnames(seu@meta.data) && "ind" %in% colnames(seu@meta.data))
  seu$donor <- seu$ind

# Step 1 — re-PCA, then Harmony-correct over the sample (stim) variable
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- RunHarmony(seu, group.by.vars = "stim",
                  reduction.use = "pca", reduction.save = "harmony", verbose = FALSE)

# Re-cluster + UMAP on the harmony embedding
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30,
               reduction.name = "umap_harmony", verbose = FALSE)

saveRDS(seu, file.path(OBJ_DIR, "ifnb_integrated.rds"))
cat("Wrote", file.path(OBJ_DIR, "ifnb_integrated.rds"), "\n")
