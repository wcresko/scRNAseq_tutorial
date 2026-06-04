#!/usr/bin/env Rscript
# Talapas analysis pipeline 07 — parallels laptop Tutorial 07 (Functional Analysis).
# Learning notebook: Exercise_Folder/Tutorial_07_FunctionalAnalysis.qmd
# Run:  sbatch --job-name=fa --time=02:00:00 --mem=32G run_rscript.sbatch 07_functional_analysis.R
# In:   ../objects/ifnb_pseudobulk_de.csv   Out: ../objects/functional/*.csv

suppressPackageStartupMessages({
  library(tidyverse); library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
})
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
OUT_DIR <- file.path(OBJ_DIR, "functional")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

de <- read_csv(file.path(OBJ_DIR, "ifnb_pseudobulk_de.csv"), show_col_types = FALSE)

# Step 2 — GO over-representation (BP) per cell type
ego_one <- function(de_ct) {
  sig <- de_ct |> filter(padj < 0.05, abs(log2FoldChange) > 1) |> pull(gene)
  if (length(sig) < 10) return(NULL)
  enrichGO(gene = sig, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
           ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = FALSE)
}

# Step 4 — GSEA on the ranked list per cell type (dedupe gene IDs first)
gsea_one <- function(de_ct) {
  rk <- de_ct |> filter(!is.na(log2FoldChange)) |>
    arrange(desc(log2FoldChange)) |> distinct(gene, .keep_all = TRUE)
  ranks <- setNames(rk$log2FoldChange, rk$gene)
  if (length(ranks) < 200) return(NULL)
  gseGO(geneList = ranks, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
        ont = "BP", pAdjustMethod = "BH", verbose = FALSE)
}

for (ct in unique(de$celltype)) {
  de_ct <- de |> filter(celltype == ct)
  safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)

  ego <- tryCatch(ego_one(de_ct), error = function(e) NULL)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0)
    write_csv(as.data.frame(ego), file.path(OUT_DIR, paste0("GO_BP_", safe_ct, ".csv")))

  gs <- tryCatch(gsea_one(de_ct), error = function(e) NULL)
  if (!is.null(gs) && nrow(as.data.frame(gs)) > 0)
    write_csv(as.data.frame(gs), file.path(OUT_DIR, paste0("GSEA_BP_", safe_ct, ".csv")))
  cat("Done:", ct, "\n")
}
