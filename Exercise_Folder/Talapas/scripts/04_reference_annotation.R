#!/usr/bin/env Rscript
# Talapas analysis pipeline 04 — parallels laptop Tutorial 04 (Reference Annotation, SingleR).
# Learning notebook: Exercise_Folder/Tutorial_04_Reference_Annotation.qmd
# Run:  sbatch --job-name=refannot run_rscript.sbatch 04_reference_annotation.R
# In:   ../objects/ifnb_annotated.rds   Out: ../objects/ifnb_annotated_final.rds
# Figures/tables (match the Mod4 notebook filenames): ../output/Mod4/Mod4_C*_*

suppressPackageStartupMessages({
  library(Seurat); library(SingleR); library(celldex)
  library(SingleCellExperiment); library(tidyverse); library(patchwork)
})
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod4") # figures/tables, named to match Tutorial_04.qmd
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
seu <- readRDS(file.path(OBJ_DIR, "ifnb_annotated.rds"))
DefaultAssay(seu) <- "RNA"

# Figure out: starting point — manual/author annotations on the UMAP (Mod4_C1)
p_setup_umap <- DimPlot(seu, group.by = "seurat_annotations", reduction = "umap",
                        label = TRUE, repel = TRUE) + NoLegend() +
  labs(title    = "Manual cell-type annotations on the unintegrated UMAP",
       subtitle = "Author-curated seurat_annotations carried from Module 3",
       x = "UMAP 1", y = "UMAP 2")
ggsave(file.path(OUT_DIR, "Mod4_C1_manual_annotations_umap.png"), p_setup_umap,
       width = 7, height = 6, dpi = 300)

# SingleR vs a celldex reference. ifnb is PBMCs, so an immune reference like the
# Human Primary Cell Atlas fits well and most cells get confident labels.
ref  <- celldex::HumanPrimaryCellAtlasData()
pred <- SingleR(test = as.SingleCellExperiment(seu), ref = ref, labels = ref$label.main)
seu$singler_label  <- pred$labels
seu$singler_pruned <- pred$pruned.labels   # NA where confidence too low
seu$singler_delta  <- pred$delta.next

# Table out: number of cells assigned to each SingleR label (Mod4_C2)
enframe(table(pred$labels), name = "singler_label", value = "n_cells") |>
  mutate(n_cells = as.integer(n_cells)) |>
  arrange(desc(n_cells)) |>
  write_csv(file.path(OUT_DIR, "Mod4_C2_singler_label_counts.csv"))

# Figure out: SingleR labels projected onto the UMAP (Mod4_C3)
p_singler_umap <- DimPlot(seu, group.by = "singler_label", reduction = "umap",
                          label = TRUE, repel = TRUE) + NoLegend() +
  labs(title    = "SingleR labels on the unintegrated UMAP",
       subtitle = "celldex HumanPrimaryCellAtlas reference (label.main)",
       x = "UMAP 1", y = "UMAP 2")
ggsave(file.path(OUT_DIR, "Mod4_C3_singler_labels_umap.png"), p_singler_umap,
       width = 7, height = 6, dpi = 300)

# Figure out: SingleR confidence diagnostics (Mod4_C4) — per-label score heatmap
# and the per-cell delta distribution (small delta = low-confidence call)
p_score_heatmap <- plotScoreHeatmap(pred)
ggsave(file.path(OUT_DIR, "Mod4_C4_singler_score_heatmap.png"), p_score_heatmap,
       width = 9, height = 7, dpi = 300)
p_delta <- plotDeltaDistribution(pred) +
  labs(title    = "SingleR delta distribution by assigned label",
       subtitle = "delta = top-label score minus median score; small delta = low-confidence call",
       x = "Assigned label", y = "Delta (score gap)")
ggsave(file.path(OUT_DIR, "Mod4_C4_singler_delta_distribution.png"), p_delta,
       width = 9, height = 5, dpi = 300)

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

# Table out: author ground-truth vs SingleR (pruned) cross-tabulation (Mod4_C7)
# (the qmd also tabulates truth vs Azimuth; this script does not run Azimuth — see header)
as.data.frame(
  table(seurat_annotations   = seu$seurat_annotations,
        singler_label_pruned = seu$singler_pruned, useNA = "ifany"),
  responseName = "n_cells") |>
  write_csv(file.path(OUT_DIR, "Mod4_C7_truth_vs_singler.csv"))

# Table out: final reconciled cell-type call counts (Mod4_C9)
enframe(table(seu$celltype_final), name = "celltype_final", value = "n_cells") |>
  mutate(n_cells = as.integer(n_cells)) |>
  arrange(desc(n_cells)) |>
  write_csv(file.path(OUT_DIR, "Mod4_C9_celltype_final_counts.csv"))

saveRDS(seu, file.path(OBJ_DIR, "ifnb_annotated_final.rds"))
message("Wrote ", file.path(OBJ_DIR, "ifnb_annotated_final.rds"))
cat("Wrote Mod4 figures/tables to", OUT_DIR, "\n")
