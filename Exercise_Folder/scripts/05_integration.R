#!/usr/bin/env Rscript
# Talapas analysis pipeline 05 — parallels laptop Tutorial 05 (Integration, Harmony).
# Learning notebook: Exercise_Folder/Tutorial_05_Integration.qmd
# Run:  sbatch --job-name=integrate run_rscript.sbatch 05_integration.R
# In:   ../data/ifnb_annotated_final.rds
# Out:  ../data/ifnb_integrated.rds  +  ../output/Mod5 figures/tables (match Tutorial_05.qmd)

suppressPackageStartupMessages({ library(Seurat); library(harmony); library(tidyverse); library(patchwork) })
set.seed(2026)
DATA_DIR <- Sys.getenv("DATA_DIR", "../data")
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod5") # figures/tables, named to match Tutorial_05.qmd
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)   # pipeline hand-off objects (.rds/.csv)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
message("[dirs] data -> ", normalizePath(DATA_DIR), "  |  figures/tables -> ", normalizePath(OUT_DIR))

# --- Figure saving --------------------------------------------------------
# save_fig() writes every figure as a .png (for viewing) AND a .svg (vector,
# editable in Illustrator / Inkscape for a manuscript). Pass the .png path; the
# .svg is written alongside with the same basename. (.svg uses the 'svglite'
# package, installed in Tutorial 00.)
save_fig <- function(filename, plot, width, height, dpi = 300, ...) {
  ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, ...)
  svg_path <- paste0(tools::file_path_sans_ext(filename), ".svg")
  tryCatch(ggplot2::ggsave(svg_path, plot, width = width, height = height, ...),
           error = function(e) message("  (could not write ", basename(svg_path),
                                       " - install 'svglite'? ", conditionMessage(e), ")"))
}
# Prefer the reference-reconciled object from step 04; fall back to step 03.
in_file <- file.path(DATA_DIR, "ifnb_annotated_final.rds")
if (!file.exists(in_file)) in_file <- file.path(DATA_DIR, "ifnb_annotated.rds")
seu <- readRDS(in_file)
DefaultAssay(seu) <- "RNA"

# ifnb's batch / sample axis is `stim` (CTRL vs IFN-beta STIM) — the two samples we
# integrate. `donor` (the 8 real `ind` patients) is carried for downstream DE / DA.
if (!"donor" %in% colnames(seu@meta.data) && "ind" %in% colnames(seu@meta.data))
  seu$donor <- seu$ind

# Step 1 — re-PCA, then Harmony-correct over the sample (stim) variable
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- RunHarmony(seu, group.by.vars = "stim",
                  reduction.use = "pca", reduction.save = "harmony", verbose = FALSE)

# Re-cluster + UMAP on the harmony embedding
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30,
               reduction.name = "umap_harmony", verbose = FALSE)

saveRDS(seu, file.path(DATA_DIR, "ifnb_integrated.rds"))
cat("Wrote", file.path(DATA_DIR, "ifnb_integrated.rds"), "\n")

# ===========================================================================
# Mod5 figures/tables — reuse the Tutorial_05.qmd filenames.
# NOTE on drift: the qmd's C6 figures compare THREE embeddings
# (unintegrated / Harmony / Seurat-anchor). This condensed pipeline computes
# Harmony only, so the saved C6 figures show the Harmony result (the by-sample
# mixing check and the by-cell-type structure check) rather than the 3-panel
# comparison. Same filenames, Harmony-only content.
# ===========================================================================

# --- Figure: sample mixing on the Harmony UMAP (Mod5_C6) -------------------
p_compare_sample <- DimPlot(seu, reduction = "umap_harmony", group.by = "stim") +
  ggtitle("Harmony — by sample") +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Sample") +
  plot_annotation(
    title    = "Sample mixing after Harmony integration",
    subtitle = "CTRL vs STIM should interleave within clusters after integration",
    caption  = "Module 5 · Two-Sample Integration (Talapas pipeline, Harmony only)")
save_fig(file.path(OUT_DIR, "Mod5_C6_compare_by_sample.png"), p_compare_sample,
       width = 7, height = 5, dpi = 300)

# --- Figure: cell-type structure on the Harmony UMAP (Mod5_C6) -------------
p_compare_cell <- DimPlot(seu, reduction = "umap_harmony",
                          group.by = "seurat_annotations", label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("Harmony — by cell type") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  plot_annotation(
    title    = "Cell-type structure after Harmony integration",
    subtitle = "Each cell type should collapse from two per-sample islands into one region",
    caption  = "Module 5 · Two-Sample Integration (Talapas pipeline, Harmony only)")
save_fig(file.path(OUT_DIR, "Mod5_C6_compare_by_celltype.png"), p_compare_cell,
       width = 7, height = 5, dpi = 300)

# --- Step 3 — Diagnose over-correction with a per-cell ISG module score -----
# Cheap call: score each cell on its IFN-beta response (scored on the RNA assay).
isg_set <- c("ISG15","IFI6","IFIT1","IFIT3","ISG20","MX1","OAS1","RSAD2","IFIT2")
DefaultAssay(seu) <- "RNA"
seu <- AddModuleScore(seu, features = list(isg_set), name = "ISG_score")

# Figure: per-cell ISG module score on the Harmony UMAP (Mod5_C7)
p_isg <- FeaturePlot(seu, features = "ISG_score1", reduction = "umap_harmony") +
  ggtitle("Harmony UMAP, ISG score (per-cell)") +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "ISG score") +
  plot_annotation(
    title    = "Per-cell IFN-beta response (ISG module score) after integration",
    subtitle = "Scored on the RNA assay; high-ISG STIM cells should co-cluster with low-ISG CTRL cells",
    caption  = "Module 5 · Two-Sample Integration")
save_fig(file.path(OUT_DIR, "Mod5_C7_isg_featureplots.png"), p_isg,
       width = 7, height = 5, dpi = 300)

# Table: mean ISG score per cell type, STIM vs CTRL (Mod5_C7)
isg_by_celltype <- seu@meta.data |>
  group_by(seurat_annotations, stim) |>
  summarise(ISG_mean = mean(ISG_score1), .groups = "drop") |>
  pivot_wider(names_from = stim, values_from = ISG_mean) |>
  arrange(desc(STIM - CTRL))
write_csv(isg_by_celltype, file.path(OUT_DIR, "Mod5_C7_isg_mean_by_celltype.csv"))

# --- Step 4 — Re-find markers on the integrated (Harmony) clusters (Mod5_C8) -
# Marker tests run on the RNA assay, not the integrated/corrected values.
# Defensive: FindAllMarkers needs a single joined RNA layer (Seurat v5 merge/split
# can leave per-sample layers); join if needed so no test is silently skipped.
if (length(Layers(seu, search = "data")) > 1) seu <- JoinLayers(seu)
Idents(seu) <- "seurat_clusters"
integrated_markers <- FindAllMarkers(seu, only.pos = TRUE,
                                     min.pct = 0.25, logfc.threshold = 0.25)
write_csv(integrated_markers, file.path(OUT_DIR, "Mod5_C8_integrated_markers.csv"))

integrated_top5 <- integrated_markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 5) |>
  ungroup()
write_csv(integrated_top5, file.path(OUT_DIR, "Mod5_C8_integrated_top5_markers.csv"))

cat("Wrote Mod5 figures/tables to", OUT_DIR, "\n")
