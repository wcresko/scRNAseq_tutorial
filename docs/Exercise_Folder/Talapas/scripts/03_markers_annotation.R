#!/usr/bin/env Rscript
# Talapas analysis pipeline 03 — parallels laptop Tutorial 03 (Markers & Annotation).
# Learning notebook: Exercise_Folder/Tutorial_03_Markers_Annotation.qmd
# Run:  sbatch --job-name=markers run_rscript.sbatch 03_markers_annotation.R
# In:   ../objects/nsclc_clustered.rds   Out: ../objects/nsclc_annotated.rds

suppressPackageStartupMessages({ library(Seurat); library(tidyverse) })
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "nsclc_clustered.rds"))

# Step 2 — Markers for every cluster (top 5 per cluster written for inspection)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
markers |> group_by(cluster) |> slice_max(order_by = avg_log2FC, n = 5) |> ungroup() |>
  write_csv(file.path(OBJ_DIR, "nsclc_top5_markers.csv"))

# Step 4 — Manual annotation by canonical markers (edit after inspecting the CSV)
canonical <- list(
  "T cells"     = c("CD3D","CD3E","TRAC"),  "B cells"     = c("MS4A1","CD79A","CD19"),
  "NK cells"    = c("NKG7","GNLY","KLRD1"), "Myeloid"     = c("LYZ","CD68","S100A8"),
  "Epithelial"  = c("EPCAM","KRT8","KRT18"),"Fibroblast"  = c("COL1A1","COL3A1","DCN"),
  "Endothelial" = c("PECAM1","VWF","CDH5"))
score_cluster <- function(cid, m) {
  cells <- WhichCells(seu, idents = cid)
  mean(rowMeans(GetAssayData(seu, layer = "data")[
    intersect(m, rownames(seu)), cells, drop = FALSE]))
}
annotation <- sapply(levels(Idents(seu)), function(cid) {
  s <- sapply(canonical, function(m) score_cluster(cid, m)); names(s)[which.max(s)]
})
seu$celltype_manual <- annotation[as.character(Idents(seu))]

# Step 6 — Save
saveRDS(seu, file.path(OBJ_DIR, "nsclc_annotated.rds"))
cat("Wrote", file.path(OBJ_DIR, "nsclc_annotated.rds"), "with",
    length(unique(seu$celltype_manual)), "cell-type labels\n")
