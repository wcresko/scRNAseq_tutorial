#!/usr/bin/env Rscript
# Talapas analysis pipeline 13 — parallels laptop Tutorial 13 (FAIR & Metadata).
# Learning notebook: Exercise_Folder/Tutorial_13_FAIR_Metadata.qmd
# Runs on the annotated NSCLC object (the laptop notebook uses the annotated ifnb).
# Run:  sbatch --job-name=fair run_rscript.sbatch 13_fair_metadata.R
# In:   ../objects/nsclc_annotated_final.rds   Out: ../objects/nsclc_for_submission.h5ad

suppressPackageStartupMessages({ library(Seurat); library(tidyverse) })
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "nsclc_annotated_final.rds"))

# Step 1 — Audit the metadata that already exists
cat("Per-cell metadata columns:\n"); print(colnames(seu@meta.data))

# Step 3 — Attach CELLxGENE-style study-level metadata (controlled-vocabulary IDs;
# edit these to match your study — look terms up at the OLS / BioPortal).
seu$organism_ontology_term_id <- "NCBITaxon:9606"   # human
seu$tissue_ontology_term_id   <- "UBERON:0002048"   # lung — set correctly for your sample
seu$assay_ontology_term_id    <- "EFO:0009922"      # 10x 3' v3
seu$suspension_type           <- "cell"
seu$is_primary_data           <- TRUE

# Step 5 — Convert Seurat -> AnnData (.h5ad). SeuratDisk's SaveH5Seurat/Convert is the
# classic bridge but can fail on some Seurat v5 objects; sceasy::convertFormat is a
# robust alternative (see the notebook).
out_h5ad <- file.path(OBJ_DIR, "nsclc_for_submission.h5ad")
if (requireNamespace("sceasy", quietly = TRUE)) {
  sceasy::convertFormat(seu, from = "seurat", to = "anndata", outFile = out_h5ad)
} else if (requireNamespace("SeuratDisk", quietly = TRUE)) {
  h5s <- file.path(OBJ_DIR, "nsclc_for_submission.h5Seurat")
  SeuratDisk::SaveH5Seurat(seu, filename = h5s, overwrite = TRUE)
  SeuratDisk::Convert(h5s, dest = "h5ad", overwrite = TRUE)
} else {
  message("Install 'sceasy' or 'SeuratDisk' to write .h5ad; see the notebook.")
}

# Step 6 — Validate against the CELLxGENE schema (shell):
#   pip install cellxgene-schema && cellxgene-schema validate ../objects/nsclc_for_submission.h5ad
cat("FAIR metadata attached; see notebook Steps 6-8 for schema validation + upload TSVs.\n")
