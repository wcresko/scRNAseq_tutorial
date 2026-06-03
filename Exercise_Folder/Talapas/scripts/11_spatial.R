#!/usr/bin/env Rscript
# Talapas analysis pipeline 11 — parallels laptop Tutorial 11 (Spatial / Visium).
# Standalone bonus track. Dataset: stxBrain (10x Visium adult mouse brain) via SeuratData.
# Learning notebook: Exercise_Folder/Tutorial_11_Spatial_Transcriptomics.qmd
# Run:  sbatch --job-name=spatial --mem=64G run_rscript.sbatch 11_spatial.R
# Out:  ../objects/brain_spatial_integrated.rds

suppressPackageStartupMessages({
  library(Seurat); library(SeuratData); library(tidyverse)
})
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
dir.create(OBJ_DIR, showWarnings = FALSE, recursive = TRUE)

# Step 1 — Load both sections (run SeuratData::InstallData("stxBrain") once first)
ant  <- LoadData("stxBrain", type = "anterior1")
post <- LoadData("stxBrain", type = "posterior1")

# Step 2 — Spatial QC metric (mouse mito prefix "^mt-")
ant[["percent.mt"]]  <- PercentageFeatureSet(ant,  pattern = "^mt-")
post[["percent.mt"]] <- PercentageFeatureSet(post, pattern = "^mt-")

# Step 3 — SCTransform normalization (recommended for Visium's wider dynamic range)
ant  <- SCTransform(ant,  assay = "Spatial", verbose = FALSE)
post <- SCTransform(post, assay = "Spatial", verbose = FALSE)

# Step 6 — Spatially variable features on the anterior section (Moran's I)
ant <- FindSpatiallyVariableFeatures(ant, assay = "SCT",
         features = VariableFeatures(ant)[1:1000], selection.method = "moransi")

# Step 7 — Merge anterior + posterior and joint-cluster (SCT integration)
ant$slice <- "anterior"; post$slice <- "posterior"
brain.list <- list(ant = ant, post = post)
features   <- SelectIntegrationFeatures(brain.list, nfeatures = 3000, verbose = FALSE)
brain.list <- PrepSCTIntegration(brain.list, anchor.features = features, verbose = FALSE)
anchors    <- FindIntegrationAnchors(brain.list, normalization.method = "SCT",
                                     anchor.features = features, verbose = FALSE)
brain      <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
DefaultAssay(brain) <- "integrated"
brain <- brain |> RunPCA(verbose = FALSE) |> FindNeighbors(dims = 1:30) |>
  FindClusters(resolution = 0.5, verbose = FALSE) |> RunUMAP(dims = 1:30)

saveRDS(brain, file.path(OBJ_DIR, "brain_spatial_integrated.rds"))
cat("Wrote", file.path(OBJ_DIR, "brain_spatial_integrated.rds"), "with", ncol(brain), "spots\n")
