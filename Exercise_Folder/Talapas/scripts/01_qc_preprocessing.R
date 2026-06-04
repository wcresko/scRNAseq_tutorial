#!/usr/bin/env Rscript
# Talapas analysis pipeline 01 — parallels laptop Tutorial 01 (QC & Preprocessing).
# Dataset: ifnb (Kang et al. 2017) — the SAME dataset as the laptop tutorials.
# Learning notebook: Exercise_Folder/Tutorial_01_QC_Preprocessing.qmd
# Run:  sbatch --job-name=qc run_rscript.sbatch 01_qc_preprocessing.R
# In:   muscData::Kang18_8vs8() (ExperimentHub cache)   Out: ../objects/ifnb_preprocessed.rds

suppressPackageStartupMessages({
  library(Seurat); library(muscData); library(SingleCellExperiment); library(tidyverse)
})
set.seed(2026)

OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")   # every step reads/writes here
dir.create(OBJ_DIR, showWarnings = FALSE, recursive = TRUE)

# Step 1 — Load ifnb from Bioconductor's ExperimentHub (cached after the first
# download). Keep singlets that have an author-assigned cell-type label, and mirror
# the metadata names the rest of the series expects:
#   - `stim`               CTRL / STIM condition
#   - `seurat_annotations` author-curated cell types
#   - `donor` (= `ind`)    the 8 real lupus donors used downstream for DE / DA
sce <- Kang18_8vs8()
sce <- sce[, sce$multiplets == "singlet" & !is.na(sce$cell)]
seu <- CreateSeuratObject(counts    = counts(sce),
                          meta.data = as.data.frame(colData(sce)),
                          project   = "ifnb", min.cells = 3, min.features = 200)
seu$stim               <- factor(toupper(seu$stim), levels = c("CTRL", "STIM"))
seu$seurat_annotations <- factor(seu$cell)
seu$donor              <- seu$ind   # real biological replicate

# Step 3 — QC metrics (SoupX needs a raw matrix and is illustrative-only on ifnb;
# scDblFinder is optional — see the notebook).
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Step 7 — Filter (PBMC-appropriate thresholds)
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Step 8 — Normalize
seu <- NormalizeData(seu)

# Step 9 — Highly variable features
seu <- FindVariableFeatures(seu, nfeatures = 2000)

# Step 10 — Scale
seu <- ScaleData(seu, vars.to.regress = c("nCount_RNA", "percent.mt"))

saveRDS(seu, file.path(OBJ_DIR, "ifnb_preprocessed.rds"))
cat("Wrote", file.path(OBJ_DIR, "ifnb_preprocessed.rds"), "with", ncol(seu), "cells\n")
