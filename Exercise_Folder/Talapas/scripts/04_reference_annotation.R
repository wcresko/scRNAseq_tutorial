#!/usr/bin/env Rscript
# Talapas analysis pipeline 04 — parallels laptop Tutorial 04 (Reference Annotation, SingleR).
# Learning notebook: Exercise_Folder/Tutorial_04_Reference_Annotation.qmd
# Run:  sbatch --job-name=refannot run_rscript.sbatch 04_reference_annotation.R
# In:   ../objects/ifnb_annotated.rds   Out: ../objects/ifnb_annotated_final.rds

suppressPackageStartupMessages({
  library(Seurat); library(SingleR); library(celldex)
  library(SingleCellExperiment); library(tidyverse)
})
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "ifnb_annotated.rds"))
DefaultAssay(seu) <- "RNA"

# SingleR vs a celldex reference. ifnb is PBMCs, so an immune reference like the
# Human Primary Cell Atlas fits well and most cells get confident labels.
ref  <- celldex::HumanPrimaryCellAtlasData()
pred <- SingleR(test = as.SingleCellExperiment(seu), ref = ref, labels = ref$label.main)
seu$singler_label  <- pred$labels
seu$singler_pruned <- pred$pruned.labels   # NA where confidence too low
seu$singler_delta  <- pred$delta.next

# Reconcile: keep the manual label where SingleR is uncertain (pruned == NA).
# Script 03 writes the manual label to `celltype_manual`.
manual_col <- intersect(c("celltype_manual","celltype","seurat_annotations"),
                        colnames(seu@meta.data))[1]
if (is.na(manual_col))
  stop("No manual-label column found (expected 'celltype_manual' from script 03).")
seu$celltype_final  <- ifelse(is.na(seu$singler_pruned),
                              as.character(seu[[manual_col]][, 1]), seu$singler_pruned)
seu$celltype_method <- ifelse(is.na(seu$singler_pruned), "manual", "singler")
print(table(manual = seu[[manual_col]][, 1], singler = seu$singler_label))

saveRDS(seu, file.path(OBJ_DIR, "ifnb_annotated_final.rds"))
message("Wrote ", file.path(OBJ_DIR, "ifnb_annotated_final.rds"))
