#!/usr/bin/env Rscript
# Talapas 09 — Functional Analysis with clusterProfiler
# Reads nsclc_pseudobulk_de.csv (from script 08) and runs:
#   - GO over-representation (BP) per cell type
#   - GSEA on a ranked list per cell type
# Writes one CSV per cell type into ../objects/functional/.

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

set.seed(2026)

OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
OUT_DIR <- file.path(OBJ_DIR, "functional")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

de <- read_csv(file.path(OBJ_DIR, "nsclc_pseudobulk_de.csv"),
               show_col_types = FALSE)

ego_one <- function(de_ct) {
  sig <- de_ct |> filter(padj < 0.05, abs(log2FoldChange) > 1) |> pull(gene)
  if (length(sig) < 10) return(NULL)
  enrichGO(gene = sig, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
           ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,
           readable = FALSE)
}

gsea_one <- function(de_ct) {
  rk <- de_ct |> filter(!is.na(log2FoldChange)) |>
    arrange(desc(log2FoldChange)) |>
    distinct(gene, .keep_all = TRUE)
  ranks <- setNames(rk$log2FoldChange, rk$gene)
  if (length(ranks) < 200) return(NULL)
  gseGO(geneList = ranks, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
        ont = "BP", pAdjustMethod = "BH", verbose = FALSE)
}

for (ct in unique(de$celltype)) {
  de_ct <- de |> filter(celltype == ct)
  safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)

  ego <- tryCatch(ego_one(de_ct), error = function(e) NULL)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write_csv(as.data.frame(ego),
              file.path(OUT_DIR, paste0("GO_BP_", safe_ct, ".csv")))
  }

  gs <- tryCatch(gsea_one(de_ct), error = function(e) NULL)
  if (!is.null(gs) && nrow(as.data.frame(gs)) > 0) {
    write_csv(as.data.frame(gs),
              file.path(OUT_DIR, paste0("GSEA_BP_", safe_ct, ".csv")))
  }
  cat("Done:", ct, "\n")
}
