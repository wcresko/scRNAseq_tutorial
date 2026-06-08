#!/usr/bin/env Rscript
# Talapas analysis pipeline 17 — parallels laptop Tutorial 17 (FAIR & Metadata).
# Learning notebook: Exercise_Folder/Tutorial_17_FAIR_Metadata.qmd
# Runs on the annotated ifnb object (the same dataset as the laptop notebook).
# Run:  sbatch --job-name=fair run_rscript.sbatch 17_fair_metadata.R
# In:   ../objects/ifnb_annotated_final.rds   Out: ../objects/ifnb_for_submission.h5ad
# Tables (match the Mod17 notebook filenames): ../output/Mod17/Mod17_C*_*

suppressPackageStartupMessages({ library(Seurat); library(tidyverse); library(patchwork) })
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod17") # tables, named to match Tutorial_17.qmd
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
seu <- readRDS(file.path(OBJ_DIR, "ifnb_annotated_final.rds"))

# Step 1 — Audit the metadata that already exists
cat("Per-cell metadata columns:\n"); print(colnames(seu@meta.data))

# Tables out: metadata column inventory and per-sample cell counts (Mod17_C2)
tibble(metadata_column = colnames(seu@meta.data)) |>
  write_csv(file.path(OUT_DIR, "Mod17_C2_metadata_columns.csv"))
enframe(table(seu$stim), name = "stim", value = "n_cells") |>
  mutate(n_cells = as.integer(n_cells)) |>
  write_csv(file.path(OUT_DIR, "Mod17_C2_cells_per_sample.csv"))

# Step 3 — Attach CELLxGENE-style study-level metadata (controlled-vocabulary IDs;
# edit these to match your study — look terms up at the OLS / BioPortal).
seu$organism_ontology_term_id <- "NCBITaxon:9606"   # human
seu$tissue_ontology_term_id   <- "UBERON:0000178"   # blood (ifnb = PBMCs)
seu$assay_ontology_term_id    <- "EFO:0009922"      # 10x 3' — ifnb is v1; set to your exact assay
seu$suspension_type           <- "cell"
seu$is_primary_data           <- TRUE

# Step 5 — Convert Seurat -> AnnData (.h5ad). SeuratDisk's SaveH5Seurat/Convert is the
# classic bridge but can fail on some Seurat v5 objects; sceasy::convertFormat is a
# robust alternative (see the notebook).
out_h5ad <- file.path(OBJ_DIR, "ifnb_for_submission.h5ad")
if (requireNamespace("sceasy", quietly = TRUE)) {
  sceasy::convertFormat(seu, from = "seurat", to = "anndata", outFile = out_h5ad)
} else if (requireNamespace("SeuratDisk", quietly = TRUE)) {
  h5s <- file.path(OBJ_DIR, "ifnb_for_submission.h5Seurat")
  SeuratDisk::SaveH5Seurat(seu, filename = h5s, overwrite = TRUE)
  SeuratDisk::Convert(h5s, dest = "h5ad", overwrite = TRUE)
} else {
  message("Install 'sceasy' or 'SeuratDisk' to write .h5ad; see the notebook.")
}

# Step 6 — Validate against the CELLxGENE schema (shell):
#   pip install cellxgene-schema && cellxgene-schema validate ../objects/ifnb_for_submission.h5ad
cat("FAIR metadata attached; see notebook Steps 6-8 for schema validation + upload TSVs.\n")
