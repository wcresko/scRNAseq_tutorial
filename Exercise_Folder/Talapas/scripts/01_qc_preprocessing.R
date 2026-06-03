#!/usr/bin/env Rscript
# Talapas analysis pipeline 01 — parallels laptop Tutorial 01 (QC & Preprocessing).
# Dataset: 20k NSCLC DTC.  Learning notebook: Exercise_Folder/Tutorial_01_QC_Preprocessing.qmd
# Run:  sbatch --job-name=qc run_rscript.sbatch 01_qc_preprocessing.R
# In:   ../data/filtered_feature_bc_matrix/   Out: ../objects/nsclc_preprocessed.rds

suppressPackageStartupMessages({ library(Seurat); library(tidyverse) })
set.seed(2026)

DATA_DIR <- Sys.getenv("DATA_DIR", "../data/filtered_feature_bc_matrix/")
OBJ_DIR  <- Sys.getenv("OBJ_DIR",  "../objects")   # every step reads/writes here
dir.create(OBJ_DIR, showWarnings = FALSE, recursive = TRUE)

# Step 1 — Load 10x counts
counts <- Read10X(data.dir = DATA_DIR)
if (is.list(counts)) counts <- counts[["Gene Expression"]]
seu <- CreateSeuratObject(counts, project = "nsclc20k", min.cells = 3, min.features = 200)

# Step 3 — QC metrics (Steps 5–6 SoupX / scDblFinder are optional; see the notebook)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Step 7 — Filter (tumor-appropriate thresholds: mt < 15)
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                            nCount_RNA < 50000 & percent.mt < 15)

# Step 8 — Normalize
seu <- NormalizeData(seu)

# Step 9 — Highly variable features
seu <- FindVariableFeatures(seu, nfeatures = 2000)

# Step 10 — Scale
seu <- ScaleData(seu, vars.to.regress = c("nCount_RNA", "percent.mt"))

saveRDS(seu, file.path(OBJ_DIR, "nsclc_preprocessed.rds"))
cat("Wrote", file.path(OBJ_DIR, "nsclc_preprocessed.rds"), "with", ncol(seu), "cells\n")
