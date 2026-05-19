#!/usr/bin/env Rscript
# Talapas 09 — Differential Abundance with miloR
# Tests for shifts in *cell-type composition* across conditions using
# neighborhood-level abundance testing on the integrated kNN graph.
# Writes nsclc_milo_da.csv.

suppressPackageStartupMessages({
  library(Seurat)
  library(miloR)
  library(SingleCellExperiment)
  library(scater)
  library(tidyverse)
})

set.seed(2026)

OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "nsclc_integrated.rds"))

# Synthetic condition (same as in script 07) if none present
if (!"condition" %in% colnames(seu@meta.data)) {
  seu$condition <- ifelse(seu$donor %in% c("d1","d2","d3"), "GroupA", "GroupB")
}

# Convert to SCE, then to Milo
sce  <- as.SingleCellExperiment(seu)
milo <- Milo(sce)

# Build kNN graph from harmony reduction
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "HARMONY")
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30,
                   refined = TRUE, reduced_dims = "HARMONY")
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)),
                   sample = "donor")

design_df <- as.data.frame(colData(milo)) |>
  distinct(donor, condition) |>
  column_to_rownames("donor")

milo <- calcNhoodDistance(milo, d = 30, reduced.dim = "HARMONY")
res  <- testNhoods(milo, design = ~ condition, design.df = design_df)

write_csv(as_tibble(res), file.path(OBJ_DIR, "nsclc_milo_da.csv"))
cat("Wrote", file.path(OBJ_DIR, "nsclc_milo_da.csv"), "\n")
