#!/usr/bin/env Rscript
# Talapas analysis pipeline 08 — parallels laptop Tutorial 08 (Differential Abundance, miloR).
# Learning notebook: Exercise_Folder/Tutorial_08_DifferentialAbundance.qmd
# Run:  sbatch --job-name=da --time=04:00:00 --mem=96G run_rscript.sbatch 08_differential_abundance.R
# In:   ../objects/ifnb_integrated.rds   Out: ../objects/ifnb_milo_da.csv
#
# ifnb has a real condition (IFN-beta STIM vs CTRL) and 8 real donors (`ind`),
# so this is genuine differential-abundance testing across neighbourhoods.

suppressPackageStartupMessages({
  library(Seurat); library(miloR); library(SingleCellExperiment)
  library(scater); library(tidyverse)
})
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "ifnb_integrated.rds"))

# Real condition + replicate. miloR needs one sample ID per replicate: donor x stim.
if (!"donor" %in% colnames(seu@meta.data)) seu$donor <- seu$ind
seu$condition <- seu$stim
seu$sample    <- paste(seu$donor, seu$stim, sep = "_")

# Step 2 — convert to Milo
sce  <- as.SingleCellExperiment(seu)
milo <- Milo(sce)

# Step 3 — build kNN graph + neighbourhoods on the harmony embedding
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "HARMONY")
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, refined = TRUE,
                   reduced_dims = "HARMONY")

# Step 4 — count cells per neighbourhood per sample
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample = "sample")

# Step 5 — design + test
design_df <- as.data.frame(colData(milo)) |>
  distinct(sample, condition) |> column_to_rownames("sample")
milo <- calcNhoodDistance(milo, d = 30, reduced.dim = "HARMONY")
res  <- testNhoods(milo, design = ~ condition, design.df = design_df)

write_csv(as_tibble(res), file.path(OBJ_DIR, "ifnb_milo_da.csv"))
cat("Wrote", file.path(OBJ_DIR, "ifnb_milo_da.csv"), "\n")
