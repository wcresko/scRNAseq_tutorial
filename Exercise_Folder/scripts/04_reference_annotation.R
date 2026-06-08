#!/usr/bin/env Rscript
# Talapas analysis pipeline 04 — parallels laptop Tutorial 04 (Reference Annotation, Azimuth).
# Learning notebook: Exercise_Folder/Tutorial_04_Reference_Annotation.qmd
# Run:  sbatch --job-name=refannot run_rscript.sbatch 04_reference_annotation.R
# In:   ../data/ifnb_annotated.rds   Out: ../data/ifnb_annotated_final.rds
# Figures/tables (match the Mod4 notebook filenames): ../output/Mod4/Mod4_C*_*
#
# NOTE: RunAzimuth installs the pbmcref.SeuratData reference (~75 MB) on first
# use. On a compute node without internet, pre-install it on a login node first
# (SeuratData::InstallData("pbmcref"), or run RunAzimuth once) so it's cached.

suppressPackageStartupMessages({
  library(Seurat); library(Azimuth)
  library(tidyverse); library(patchwork)
})
set.seed(2026)
DATA_DIR <- Sys.getenv("DATA_DIR", "../data")
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod4") # figures/tables, named to match Tutorial_04.qmd
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)   # pipeline hand-off objects (.rds/.csv)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("[dirs] data -> ", normalizePath(DATA_DIR), "  |  figures/tables -> ", normalizePath(OUT_DIR))
seu <- readRDS(file.path(DATA_DIR, "ifnb_annotated.rds"))
DefaultAssay(seu) <- "RNA"

# Figure out: starting point — manual/author annotations on the UMAP (Mod4_C1)
p_setup_umap <- DimPlot(seu, group.by = "seurat_annotations", reduction = "umap",
                        label = TRUE, repel = TRUE) + NoLegend() +
  labs(title    = "Manual cell-type annotations on the unintegrated UMAP",
       subtitle = "Author-curated seurat_annotations carried from Module 3",
       x = "UMAP 1", y = "UMAP 2")
ggsave(file.path(OUT_DIR, "Mod4_C1_manual_annotations_umap.png"), p_setup_umap,
       width = 7, height = 6, dpi = 300)

# Azimuth anchor-based reference mapping. ifnb is PBMCs, so the pbmcref reference
# fits well and most cells project cleanly (high mapping score).
seu <- RunAzimuth(seu, reference = "pbmcref")

# Table out: per-cell Azimuth predictions (l1/l2/l3 + scores) (Mod4_C2)
seu@meta.data[, grep("predicted", colnames(seu@meta.data)), drop = FALSE] |>
  rownames_to_column("cell") |>
  write_csv(file.path(OUT_DIR, "Mod4_C2_azimuth_predictions.csv"))

# Figure out: Azimuth level-2 labels projected onto the UMAP (Mod4_C2)
p_azimuth_umap <- DimPlot(seu, group.by = "predicted.celltype.l2", reduction = "umap",
                          label = TRUE, repel = TRUE) + NoLegend() +
  labs(title    = "Azimuth predicted cell types (level 2) on the unintegrated UMAP",
       subtitle = "Anchor-based mapping against the Azimuth pbmcref reference",
       x = "UMAP 1", y = "UMAP 2")
ggsave(file.path(OUT_DIR, "Mod4_C2_azimuth_labels_umap.png"), p_azimuth_umap,
       width = 7, height = 6, dpi = 300)

# Figure out: Azimuth mapping score per cell (Mod4_C3) — low = off-manifold
p_mapping <- FeaturePlot(seu, features = "mapping.score", reduction = "umap") +
  labs(title    = "Azimuth mapping score per cell",
       subtitle = "High (~0.7-1.0) = projects cleanly into the reference; low = off-manifold",
       x = "UMAP 1", y = "UMAP 2")
ggsave(file.path(OUT_DIR, "Mod4_C3_azimuth_mapping_score.png"), p_mapping,
       width = 7, height = 6, dpi = 300)

# Reconcile: keep the manual/ground-truth label where Azimuth agrees, else flag.
# Script 03 writes the manual label to `celltype_manual`.
manual_col <- intersect(c("celltype_manual","celltype","seurat_annotations"),
                        colnames(seu@meta.data))[1]
if (is.na(manual_col))
  stop("No manual-label column found (expected 'celltype_manual' from script 03).")
agree <- as.character(seu[[manual_col]][, 1]) == as.character(seu$predicted.celltype.l2)
seu$celltype_final  <- ifelse(agree, as.character(seu[[manual_col]][, 1]),
                              paste0("REVIEW:", as.character(seu[[manual_col]][, 1])))
seu$celltype_method <- ifelse(agree, "consensus", "manual")
print(table(manual = seu[[manual_col]][, 1], azimuth = seu$predicted.celltype.l2))

# Table out: author ground-truth vs Azimuth cross-tabulation (Mod4_C4)
as.data.frame(
  table(seurat_annotations    = seu$seurat_annotations,
        predicted_celltype_l2 = seu$predicted.celltype.l2, useNA = "ifany"),
  responseName = "n_cells") |>
  write_csv(file.path(OUT_DIR, "Mod4_C4_truth_vs_azimuth.csv"))

# Table out: final reconciled cell-type call counts (Mod4_C6)
enframe(table(seu$celltype_final), name = "celltype_final", value = "n_cells") |>
  mutate(n_cells = as.integer(n_cells)) |>
  arrange(desc(n_cells)) |>
  write_csv(file.path(OUT_DIR, "Mod4_C6_celltype_final_counts.csv"))

saveRDS(seu, file.path(DATA_DIR, "ifnb_annotated_final.rds"))
message("Wrote ", file.path(DATA_DIR, "ifnb_annotated_final.rds"))
cat("Wrote Mod4 figures/tables to", OUT_DIR, "\n")
