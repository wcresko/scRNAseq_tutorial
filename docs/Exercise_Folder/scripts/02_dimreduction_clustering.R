#!/usr/bin/env Rscript
# Talapas analysis pipeline 02 — parallels laptop Tutorial 02 (DimReduction & Clustering).
# Learning notebook: Exercise_Folder/Tutorial_02_DimReduction_Clustering.qmd
# Run:  sbatch --job-name=cluster run_rscript.sbatch 02_dimreduction_clustering.R
# In:   ../data/ifnb_preprocessed.rds   Out: ../data/ifnb_clustered.rds
# Figures/tables (match the Mod2 notebook filenames): ../output/Mod2/Mod2_C*_*

suppressPackageStartupMessages({ library(Seurat); library(tidyverse); library(patchwork) })
set.seed(2026)
DATA_DIR <- Sys.getenv("DATA_DIR", "../data")
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod2") # figures/tables, named to match Tutorial_02.qmd
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)   # pipeline hand-off objects (.rds/.csv)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("[dirs] data -> ", normalizePath(DATA_DIR), "  |  figures/tables -> ", normalizePath(OUT_DIR))

# --- Figure saving --------------------------------------------------------
# save_fig(): write each figure as .png + .svg (vector). See Tutorial_02_DimReduction_Clustering.qmd.
save_fig <- function(filename, plot, width, height, dpi = 300, ...) {
  ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, ...)
  svg_path <- paste0(tools::file_path_sans_ext(filename), ".svg")
  tryCatch(ggplot2::ggsave(svg_path, plot, width = width, height = height, ...),
           error = function(e) message("  (could not write ", basename(svg_path),
                                       " - install 'svglite'? ", conditionMessage(e), ")"))
}
seu <- readRDS(file.path(DATA_DIR, "ifnb_preprocessed.rds"))

# Step 1 — PCA
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

# Table out: top +/- gene loadings for the first five PCs (Mod2_C3)
pca_loadings <- Loadings(seu[["pca"]])[, 1:5]
purrr::map_dfr(1:5, function(d) {
  ord <- order(pca_loadings[, d])
  tibble(PC        = paste0("PC", d),
         direction = c(rep("negative", 5), rep("positive", 5)),
         gene      = c(rownames(pca_loadings)[head(ord, 5)],
                       rownames(pca_loadings)[tail(ord, 5)]),
         loading   = c(pca_loadings[head(ord, 5), d],
                       pca_loadings[tail(ord, 5), d]))
}) |>
  write_csv(file.path(OUT_DIR, "Mod2_C3_pca_loadings.csv"))

# Figure out: PC1 gene-loading heatmap (Mod2_C3) — expect ISGs to dominate PC1
p_pca_heatmap <- DimHeatmap(seu, dims = 1, cells = 500, balanced = TRUE, fast = FALSE) +
  ggtitle("PC1 loadings — top genes across 500 cells")
p_pca_heatmap <- wrap_elements(p_pca_heatmap) + plot_annotation(
  title    = "PCA diagnostic — PC1 gene-loading heatmap",
  subtitle = "Cells ordered by PC1 score; expect interferon-stimulated genes to dominate",
  caption  = "Module 2 · Dimensionality Reduction & Clustering (Talapas pipeline)")
save_fig(file.path(OUT_DIR, "Mod2_C3_pca_heatmap.png"), p_pca_heatmap, width = 7, height = 6, dpi = 300)

# Figure out: elbow plot — variance per PC (Mod2_C4)
p_elbow <- ElbowPlot(seu, ndims = 40) +
  labs(title    = "PCA elbow plot — variance captured per principal component",
       subtitle = "The elbow (~PC15–20) marks where additional PCs stop adding signal",
       x = "Principal component", y = "Standard deviation")
save_fig(file.path(OUT_DIR, "Mod2_C4_elbow.png"), p_elbow, width = 7, height = 5, dpi = 300)

# Step 2 — Choose PCs (inspect ElbowPlot interactively). Kept in sync with the
# laptop notebook (Tutorial_02): the elbow on ifnb flattens by ~PC 20.
N_PCS <- 20

# Step 3 — Neighbor graph
seu <- FindNeighbors(seu, dims = 1:N_PCS, verbose = FALSE)

# Step 4 — Cluster at several resolutions (Louvain default; see Tutorial_02.qmd for Leiden)
for (res in c(0.1, 0.3, 0.5, 0.7, 1.0))
  seu <- FindClusters(seu, resolution = res, verbose = FALSE)
Idents(seu) <- "RNA_snn_res.0.5"

# Step 5 — UMAP
seu <- RunUMAP(seu, dims = 1:N_PCS, verbose = FALSE)

# Figure out: UMAP coloured by graph cluster, res = 0.5 (Mod2_C8)
p_umap_clusters <- DimPlot(seu, reduction = "umap", label = TRUE) +
  labs(title    = "Unintegrated UMAP — graph clusters (res = 0.5)",
       subtitle = "Clusters computed on the joint CTRL + STIM matrix without integration",
       x = "UMAP 1", y = "UMAP 2", colour = "Cluster")
save_fig(file.path(OUT_DIR, "Mod2_C8_umap_clusters.png"), p_umap_clusters, width = 7, height = 6, dpi = 300)

# Step 6 — Diagnose the batch effect: colour the SAME UMAP by sample and by cell type
p_by_sample <- DimPlot(seu, reduction = "umap", group.by = "stim") +
  labs(title    = "Unintegrated UMAP — coloured by sample (stim)",
       subtitle = "CTRL and STIM cells form largely parallel islands",
       x = "UMAP 1", y = "UMAP 2", colour = "Sample")
p_by_celltype <- DimPlot(seu, reduction = "umap", group.by = "seurat_annotations",
                         label = TRUE, repel = TRUE) + NoLegend() +
  labs(title    = "Unintegrated UMAP — coloured by author cell-type",
       subtitle = "Most cell types split into two sister clusters, tracking stim",
       x = "UMAP 1", y = "UMAP 2")
p_batch <- (p_by_sample | p_by_celltype) + plot_annotation(
  title    = "Batch-effect diagnostic — sample vs cell-type on the same embedding",
  subtitle = "If cell types co-clustered, the two panels would overlay; here they don't",
  caption  = "Module 2 · Dimensionality Reduction & Clustering (Talapas pipeline)")
save_fig(file.path(OUT_DIR, "Mod2_C9_umap_by_sample.png"),   p_by_sample,   width = 7,  height = 6, dpi = 300)
save_fig(file.path(OUT_DIR, "Mod2_C9_umap_by_celltype.png"), p_by_celltype, width = 7,  height = 6, dpi = 300)
save_fig(file.path(OUT_DIR, "Mod2_C9_umap_batch_diagnostic.png"), p_batch,  width = 13, height = 6, dpi = 300)

# Tables out: cluster composition by sample and by cell type (Mod2_C10)
as.data.frame(round(prop.table(
  table(cluster = seu$RNA_snn_res.0.5, sample = seu$stim), margin = 1), 2)) |>
  write_csv(file.path(OUT_DIR, "Mod2_C10_cluster_by_sample.csv"))
as.data.frame(
  table(cluster = seu$RNA_snn_res.0.5, celltype = seu$seurat_annotations)) |>
  write_csv(file.path(OUT_DIR, "Mod2_C10_cluster_by_celltype.csv"))

saveRDS(seu, file.path(DATA_DIR, "ifnb_clustered.rds"))
cat("Wrote", file.path(DATA_DIR, "ifnb_clustered.rds"), "with",
    length(levels(Idents(seu))), "clusters\n")
cat("Wrote Mod2 figures/tables to", OUT_DIR, "\n")
