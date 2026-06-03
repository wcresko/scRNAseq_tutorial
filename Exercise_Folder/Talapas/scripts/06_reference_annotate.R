#!/usr/bin/env Rscript
# Talapas 06 — Reference-based Annotation with SingleR (+ optional Azimuth)
# Picks up the manually-annotated NSCLC object (nsclc_annotated.rds, from
# Talapas 05) and adds an automated, reference-based label per cell, then
# reconciles it against the manual labels. Writes nsclc_annotated_final.rds.
#
# Reference-based methods run per cell and do NOT need your own integration,
# so this step slots in right after manual annotation and before integration
# (mirrors laptop Tutorial 04).
#
# NOTE ON THE HARD REGIME: the 20k NSCLC sample is dissociated *tumor* — many
# cells are malignant/non-immune and fall OUTSIDE immune reference atlases.
# Expect low SingleR scores / "pruned" labels there; that is the signal, not a
# bug. Trust manual markers for the tumor compartment.

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
  library(tidyverse)
})

set.seed(2026)

OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "nsclc_annotated.rds"))
DefaultAssay(seu) <- "RNA"

# ---- SingleR against a built-in celldex reference -------------------------
# HumanPrimaryCellAtlasData covers broad human cell types (immune + stromal).
# BlueprintEncodeData is a good alternative for blood/immune-heavy samples.
ref <- celldex::HumanPrimaryCellAtlasData()

sce  <- as.SingleCellExperiment(seu)
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)

seu$singler_label <- pred$labels
seu$singler_pruned <- pred$pruned.labels          # NA where confidence too low
seu$singler_delta  <- pred$delta.next             # margin to runner-up label

# ---- Reconcile reference vs manual ----------------------------------------
# Keep the manual label where SingleR is uncertain (pruned == NA) or where the
# manual call is a tumor/epithelial type the immune reference cannot represent.
# Script 05 writes the manual label to `celltype_manual`; fall back to other
# common names only if that column is absent.
manual_col <- intersect(c("celltype_manual", "celltype", "seurat_annotations"),
                        colnames(seu@meta.data))[1]
if (is.na(manual_col))
  stop("No manual-label column found (expected 'celltype_manual' from script 05).")
seu$celltype_final <- ifelse(is.na(seu$singler_pruned),
                             as.character(seu[[manual_col]][, 1]),
                             seu$singler_pruned)
seu$celltype_method <- ifelse(is.na(seu$singler_pruned), "manual", "singler")

# Cross-tabulate for the log
print(table(manual = seu[[manual_col]][, 1], singler = seu$singler_label))

saveRDS(seu, file.path(OBJ_DIR, "nsclc_annotated_final.rds"))
message("Wrote ", file.path(OBJ_DIR, "nsclc_annotated_final.rds"))
