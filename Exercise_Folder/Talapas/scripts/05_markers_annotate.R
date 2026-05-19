#!/usr/bin/env Rscript
# Talapas 05 — Markers & Cell Type Annotation
# Reads nsclc_clustered.rds, finds cluster markers, applies manual labels
# from a curated PBMC/NSCLC marker set. Writes nsclc_annotated.rds.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

set.seed(2026)

OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "nsclc_clustered.rds"))

# Find markers for all clusters
markers <- FindAllMarkers(seu,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25,
                          verbose = FALSE)

# Top 5 markers per cluster (for inspection)
top5 <- markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = 5) |>
  ungroup()
write_csv(top5, file.path(OBJ_DIR, "nsclc_top5_markers.csv"))

# Lightweight manual annotation by canonical markers.
# Replace with your own curated mapping after inspecting top5_markers.csv.
canonical <- list(
  "T cells"      = c("CD3D", "CD3E", "TRAC"),
  "B cells"      = c("MS4A1", "CD79A", "CD19"),
  "NK cells"     = c("NKG7", "GNLY", "KLRD1"),
  "Myeloid"      = c("LYZ", "CD68", "S100A8"),
  "Epithelial"   = c("EPCAM", "KRT8", "KRT18"),
  "Fibroblast"   = c("COL1A1", "COL3A1", "DCN"),
  "Endothelial"  = c("PECAM1", "VWF", "CDH5")
)

score_cluster <- function(seu, cluster_id, markers) {
  cells <- WhichCells(seu, idents = cluster_id)
  mean(rowMeans(GetAssayData(seu, layer = "data")[
    intersect(markers, rownames(seu)), cells, drop = FALSE]))
}

annotation <- sapply(levels(Idents(seu)), function(cid) {
  scores <- sapply(canonical, function(m) score_cluster(seu, cid, m))
  names(scores)[which.max(scores)]
})
seu$celltype_manual <- annotation[as.character(Idents(seu))]

saveRDS(seu, file.path(OBJ_DIR, "nsclc_annotated.rds"))
cat("Wrote", file.path(OBJ_DIR, "nsclc_annotated.rds"),
    "with", length(unique(seu$celltype_manual)), "cell-type labels\n")
