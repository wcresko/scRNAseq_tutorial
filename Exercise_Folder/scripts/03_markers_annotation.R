#!/usr/bin/env Rscript
# Talapas analysis pipeline 03 — parallels laptop Tutorial 03 (Markers & Annotation).
# Learning notebook: Exercise_Folder/Tutorial_03_Markers_Annotation.qmd
# Run:  sbatch --job-name=markers run_rscript.sbatch 03_markers_annotation.R
# In:   ../data/ifnb_clustered.rds   Out: ../data/ifnb_annotated.rds
# Figures/tables (match the Mod3 notebook filenames): ../output/Mod3/Mod3_C*_*

suppressPackageStartupMessages({ library(Seurat); library(tidyverse); library(patchwork) })
set.seed(2026)
DATA_DIR <- Sys.getenv("DATA_DIR", "../data")
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod3") # figures/tables, named to match Tutorial_03.qmd
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)   # pipeline hand-off objects (.rds/.csv)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
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
seu <- readRDS(file.path(DATA_DIR, "ifnb_clustered.rds"))

# Step 2 — Markers for every cluster (top 5 per cluster written for inspection)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
top5 <- markers |> group_by(cluster) |> slice_max(order_by = avg_log2FC, n = 5) |> ungroup()
top5 |> write_csv(file.path(DATA_DIR, "ifnb_top5_markers.csv"))   # pipeline copy (read downstream)

# Tables out: all positive markers + the per-cluster top 5 (Mod3_C4)
markers |> write_csv(file.path(OUT_DIR, "Mod3_C4_all_markers.csv"))
top5    |> write_csv(file.path(OUT_DIR, "Mod3_C4_top5_markers_per_cluster.csv"))

# Step 3 — Visualise canonical PBMC markers (object already carries the umap reduction)
canonical_vln <- c("CD14", "LYZ", "MS4A1", "CD79A", "CD8A", "GNLY")
p_vln <- VlnPlot(seu, features = canonical_vln, pt.size = 0, ncol = 3) &
  xlab("Cluster (RNA_snn_res.0.5)") & ylab("Log-normalized expression") &
  theme(plot.title = element_text(size = 11))
p_vln <- p_vln + plot_annotation(
  title    = "Canonical PBMC marker expression by cluster",
  subtitle = "Monocyte (CD14/LYZ), B (MS4A1/CD79A), CD8 T (CD8A) and NK (GNLY) markers",
  caption  = "Module 3 · Markers & Annotation (Talapas pipeline)")
save_fig(file.path(OUT_DIR, "Mod3_C5_marker_violins.png"), p_vln, width = 10, height = 6, dpi = 300)

p_feature <- FeaturePlot(seu, features = c("CD14", "MS4A1", "CD8A", "GNLY"), ncol = 2) &
  xlab("UMAP 1") & ylab("UMAP 2") & theme(plot.title = element_text(size = 11))
p_feature <- p_feature + plot_annotation(
  title    = "Marker expression overlaid on the UMAP",
  subtitle = "CD14 (monocytes), MS4A1 (B), CD8A (CD8 T), GNLY (NK)",
  caption  = "Module 3 · Markers & Annotation (Talapas pipeline)")
save_fig(file.path(OUT_DIR, "Mod3_C6_marker_featureplots.png"), p_feature, width = 9, height = 8, dpi = 300)

# Heatmap: re-scale the top-marker genes so non-HVG markers also appear (Tutorial 01 scaled HVGs only)
top10 <- markers |> group_by(cluster) |> slice_max(avg_log2FC, n = 10) |> ungroup()
seu <- ScaleData(seu, features = unique(top10$gene), verbose = FALSE)
p_heatmap <- DoHeatmap(seu, features = top10$gene) + NoLegend() +
  labs(title    = "Top-10 markers per cluster — scaled expression heatmap",
       subtitle = "Each column is a cell, grouped by cluster; rows are marker genes",
       caption  = "Module 3 · Markers & Annotation (Talapas pipeline)")
save_fig(file.path(OUT_DIR, "Mod3_C7_marker_heatmap.png"), p_heatmap, width = 12, height = 10, dpi = 300)

canonical_panel <- c("CD14", "LYZ", "FCGR3A", "MS4A7", "FCER1A", "CST3",
                     "IL7R", "CCR7", "S100A4", "CD8A", "GZMK", "GNLY",
                     "NKG7", "MS4A1", "CD79A", "PPBP")
p_dot <- DotPlot(seu, features = canonical_panel) + RotatedAxis() +
  labs(title    = "Canonical PBMC marker panel across clusters",
       subtitle = "Dot size = % of cells expressing; colour = scaled average expression",
       x = "Marker gene", y = "Cluster (RNA_snn_res.0.5)",
       colour = "Avg. expression", size = "Percent expressed")
save_fig(file.path(OUT_DIR, "Mod3_C8_marker_dotplot.png"), p_dot, width = 10, height = 6, dpi = 300)

# Step 4 — Manual annotation by canonical PBMC markers (edit after inspecting the CSV)
canonical <- list(
  "CD14+ Mono" = c("CD14","LYZ","S100A8"),    "CD16+ Mono" = c("FCGR3A","MS4A7"),
  "B"          = c("MS4A1","CD79A","CD79B"),  "CD4 T"      = c("IL7R","CCR7","CD3D"),
  "CD8 T"      = c("CD8A","CD8B","GZMK"),      "NK"         = c("GNLY","NKG7","KLRD1"),
  "DC"         = c("FCER1A","CST3"),           "pDC"        = c("GZMB","SERPINF1"),
  "Mk"         = c("PPBP","PF4"))
score_cluster <- function(cid, m) {
  cells <- WhichCells(seu, idents = cid)
  mean(rowMeans(GetAssayData(seu, layer = "data")[
    intersect(m, rownames(seu)), cells, drop = FALSE]))
}
annotation <- sapply(levels(Idents(seu)), function(cid) {
  s <- sapply(canonical, function(m) score_cluster(cid, m)); names(s)[which.max(s)]
})
seu$celltype_manual <- annotation[as.character(Idents(seu))]

# Step 5 — Validate manual labels against the author ground truth (Mod3_C10)
p_manual <- DimPlot(seu, reduction = "umap", group.by = "celltype_manual",
                    label = TRUE, repel = TRUE) + NoLegend() +
  labs(title = "Manual annotation (yours)", x = "UMAP 1", y = "UMAP 2")
p_truth  <- DimPlot(seu, reduction = "umap", group.by = "seurat_annotations",
                    label = TRUE, repel = TRUE) + NoLegend() +
  labs(title = "Author labels (ground truth)", x = "UMAP 1", y = "UMAP 2")
p_validate <- (p_manual | p_truth) + plot_annotation(
  title    = "Manual annotation vs author ground truth on the same embedding",
  subtitle = "Cells should land in matching positions if your manual labels agree with the authors'",
  caption  = "Module 3 · Markers & Annotation (Talapas pipeline)")
save_fig(file.path(OUT_DIR, "Mod3_C10_umap_manual.png"),          p_manual,   width = 7,  height = 6, dpi = 300)
save_fig(file.path(OUT_DIR, "Mod3_C10_umap_ground_truth.png"),    p_truth,    width = 7,  height = 6, dpi = 300)
save_fig(file.path(OUT_DIR, "Mod3_C10_umap_manual_vs_truth.png"), p_validate, width = 13, height = 6, dpi = 300)

# Table out: cross-tab of manual labels vs ground-truth labels (Mod3_C10)
as.data.frame(table(manual = seu$celltype_manual, truth = seu$seurat_annotations)) |>
  write_csv(file.path(OUT_DIR, "Mod3_C10_manual_vs_truth.csv"))

# Step 6 — Save
saveRDS(seu, file.path(DATA_DIR, "ifnb_annotated.rds"))
cat("Wrote", file.path(DATA_DIR, "ifnb_annotated.rds"), "with",
    length(unique(seu$celltype_manual)), "cell-type labels\n")
cat("Wrote Mod3 figures/tables to", OUT_DIR, "\n")
