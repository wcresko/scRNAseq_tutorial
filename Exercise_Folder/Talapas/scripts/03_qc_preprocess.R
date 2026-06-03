#!/usr/bin/env Rscript
# Talapas 03 — QC & Preprocessing on 20k NSCLC DTC
# Reads filtered_feature_bc_matrix/, applies standard Seurat QC, writes a
# preprocessed .rds for the next step.
#
# Run via:
#   sbatch --account=<PIRG> --partition=compute --time=04:00:00 \
#          --cpus-per-task=8 --mem=64G --wrap "Rscript 03_qc_preprocess.R"

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

set.seed(2026)

DATA_DIR <- Sys.getenv("DATA_DIR", "../data/filtered_feature_bc_matrix/")
OBJ_DIR  <- Sys.getenv("OBJ_DIR",  "../objects")   # every step reads/writes here
dir.create(OBJ_DIR, showWarnings = FALSE, recursive = TRUE)

# Load 10x counts
counts <- Read10X(data.dir = DATA_DIR)
if (is.list(counts)) counts <- counts[["Gene Expression"]]

seu <- CreateSeuratObject(counts, project = "nsclc20k",
                          min.cells = 3, min.features = 200)

# Mitochondrial percent
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# QC filter
seu <- subset(seu,
              subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                       nCount_RNA   < 50000 &
                       percent.mt   < 15)

# Normalize, find HVGs, scale
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu <- ScaleData(seu, vars.to.regress = c("nCount_RNA", "percent.mt"))

saveRDS(seu, file.path(OBJ_DIR, "nsclc_preprocessed.rds"))
cat("Wrote", file.path(OBJ_DIR, "nsclc_preprocessed.rds"),
    "with", ncol(seu), "cells\n")
