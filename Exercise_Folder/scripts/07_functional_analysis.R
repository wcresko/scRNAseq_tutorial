#!/usr/bin/env Rscript
# Talapas analysis pipeline 07 — parallels laptop Tutorial 07 (Functional Analysis).
# Learning notebook: Exercise_Folder/Tutorial_07_FunctionalAnalysis.qmd
# Run:  sbatch --job-name=fa --time=02:00:00 --mem=32G run_rscript.sbatch 07_functional_analysis.R
# In:   ../data/ifnb_pseudobulk_de.csv
# Out:  ../data/functional/*.csv  +  ../output/Mod7 figures/tables (match Tutorial_07.qmd)

suppressPackageStartupMessages({
  library(tidyverse); library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(ReactomePA); library(patchwork)
})
set.seed(2026)
DATA_DIR <- Sys.getenv("DATA_DIR", "../data")
FUN_DIR <- file.path(DATA_DIR, "functional")          # pipeline per-cell-type result CSVs
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod7")    # figures/tables, named to match Tutorial_07.qmd
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)   # pipeline hand-off objects (.rds/.csv)
dir.create(FUN_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("[dirs] data -> ", normalizePath(DATA_DIR), "  |  figures/tables -> ", normalizePath(OUT_DIR))

de <- read_csv(file.path(DATA_DIR, "ifnb_pseudobulk_de.csv"), show_col_types = FALSE)

# Step 1 — map gene SYMBOL -> ENTREZID (cheap bitr call; enrichment tools key on
# ENTREZ). A 5-10% unmapped rate (deprecated symbols, lncRNAs, pseudogenes) is normal.
de_mapped <- de |>
  mutate(symbol = gene) |>
  inner_join(
    bitr(unique(de$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    by = c("symbol" = "SYMBOL"))
cat("Mapped:   ", n_distinct(de_mapped$symbol), "\n")
cat("Unmapped: ", n_distinct(de$gene) - n_distinct(de_mapped$symbol), "\n")

# Step 2 — GO over-representation (BP) per cell type. The universe is the set of
# genes actually TESTED for that cell type (not the whole genome).
run_enrichGO <- function(ct) {
  sub <- de_mapped |> filter(celltype == ct)
  sig <- sub |> filter(padj < 0.05, abs(log2FoldChange) > 0.5) |> pull(ENTREZID)
  uni <- sub |> pull(ENTREZID)
  if (length(sig) < 10) return(NULL)
  ego <- enrichGO(gene = sig, universe = uni, OrgDb = org.Hs.eg.db,
                  ont = "BP", keyType = "ENTREZID", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.10, readable = TRUE)
  if (!is.null(ego)) ego@result$celltype <- ct
  ego
}
ego_list <- map(unique(de_mapped$celltype), run_enrichGO) |>
  set_names(unique(de_mapped$celltype)) |> compact()

# Table: top-5 GO BP terms per cell type (Mod7_C3)
ego_top5 <- map_dfr(names(ego_list),
  \(ct) head(ego_list[[ct]]@result, 5) |>
        mutate(celltype = ct) |>
        select(celltype, ID, Description, Count, p.adjust))
write_csv(ego_top5, file.path(OUT_DIR, "Mod7_C3_enrichgo_top5.csv"))

# Step 3 — cross-cell-type comparison with compareCluster (Mod7_C4)
sig_by_ct <- de_mapped |>
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) |>
  split(~ celltype) |>
  map(~ unique(.x$ENTREZID))
sig_by_ct <- compact(sig_by_ct[lengths(sig_by_ct) >= 10])

cc <- tryCatch(
  compareCluster(geneClusters = sig_by_ct, fun = "enrichGO", OrgDb = org.Hs.eg.db,
                 ont = "BP", keyType = "ENTREZID", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05),
  error = function(e) NULL)

# Figure: GO BP over-representation per cell type as a dot plot (Mod7_C4)
if (!is.null(cc)) {
  p_cc <- dotplot(cc, showCategory = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title    = "GO BP enrichment per cell type — STIM vs CTRL",
         subtitle = "Top 8 over-represented Biological Process terms per cell type",
         x        = "Cell type",
         y        = "GO Biological Process term",
         colour   = "Adjusted p-value",
         size     = "Gene ratio")
  ggsave(file.path(OUT_DIR, "Mod7_C4_comparecluster_dotplot.png"), p_cc,
         width = 10, height = 8, dpi = 300)
}

# Step 4 — GSEA on the ranked log2FC list per cell type (dedupe ENTREZIDs first
# to avoid the rank-tie warning).
run_gseGO <- function(ct) {
  sub <- de_mapped |>
    filter(celltype == ct, !is.na(ENTREZID)) |>
    arrange(desc(log2FoldChange)) |> distinct(ENTREZID, .keep_all = TRUE)
  rk <- setNames(sub$log2FoldChange, sub$ENTREZID)
  if (length(rk) < 200) return(NULL)
  gse <- gseGO(geneList = rk, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID",
               pAdjustMethod = "BH", pvalueCutoff = 0.05, eps = 0, verbose = FALSE)
  if (!is.null(gse)) gse@result$celltype <- ct
  gse
}
gse_list <- map(unique(de_mapped$celltype), run_gseGO) |>
  set_names(unique(de_mapped$celltype)) |> compact()

# Table: top-5 GSEA BP terms per cell type, ranked by NES (Mod7_C5)
gse_top5 <- map_dfr(names(gse_list),
  \(ct) gse_list[[ct]]@result |>
        arrange(desc(NES)) |> head(5) |>
        mutate(celltype = ct) |>
        select(celltype, Description, NES, p.adjust))
write_csv(gse_top5, file.path(OUT_DIR, "Mod7_C5_gsego_top5.csv"))

# Figure: GSEA running-enrichment plot for the top 3 BP terms in CD14 monocytes (Mod7_C6)
ct_pick <- "CD14 Mono"
if (!is.null(gse_list[[ct_pick]]) && nrow(gse_list[[ct_pick]]@result) > 0) {
  p_gsea <- gseaplot2(gse_list[[ct_pick]],
                      geneSetID = head(gse_list[[ct_pick]]@result$ID, 3),
                      title     = paste0(ct_pick, " — top 3 BP enrichments"))
  ggsave(file.path(OUT_DIR, "Mod7_C6_gsea_running_enrichment.png"), p_gsea,
         width = 8, height = 6, dpi = 300)
}

# Step 5 — Reactome pathway over-representation per cell type (cleaner than GO BP
# for immune signaling).
reactome_list <- map(unique(de_mapped$celltype), function(ct) {
  sub <- de_mapped |> filter(celltype == ct)
  sig <- sub |> filter(padj < 0.05, abs(log2FoldChange) > 0.5) |> pull(ENTREZID)
  if (length(sig) < 10) return(NULL)
  enrichPathway(gene = sig, organism = "human", pAdjustMethod = "BH",
                pvalueCutoff = 0.05, readable = TRUE)
}) |> set_names(unique(de_mapped$celltype)) |> compact()

# Table: top-5 Reactome pathways per cell type (Mod7_C7)
reactome_top5 <- map_dfr(names(reactome_list),
  \(ct) head(reactome_list[[ct]]@result, 5) |>
        mutate(celltype = ct) |>
        select(celltype, Description, p.adjust))
write_csv(reactome_top5, file.path(OUT_DIR, "Mod7_C7_reactome_top5.csv"))

# Step 6 — Save the full per-cell-type enrichment tables (Mod7_C8). Also keep the
# pipeline's per-cell-type CSV copies in DATA_DIR/functional for downstream use.
ego_long <- map_dfr(names(ego_list), \(ct) ego_list[[ct]]@result |> mutate(celltype = ct))
write_csv(ego_long, file.path(OUT_DIR, "Mod7_C8_enrichgo_bp.csv"))

reactome_long <- map_dfr(names(reactome_list),
                         \(ct) reactome_list[[ct]]@result |> mutate(celltype = ct))
write_csv(reactome_long, file.path(OUT_DIR, "Mod7_C8_reactome.csv"))

gse_long <- map_dfr(names(gse_list), \(ct) gse_list[[ct]]@result |> mutate(celltype = ct))
write_csv(gse_long, file.path(OUT_DIR, "Mod7_C8_gsego_bp.csv"))

# Pipeline per-cell-type result CSVs (kept in DATA_DIR/functional)
for (ct in names(ego_list)) {
  safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)
  write_csv(ego_list[[ct]]@result, file.path(FUN_DIR, paste0("GO_BP_", safe_ct, ".csv")))
}
for (ct in names(gse_list)) {
  safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)
  write_csv(gse_list[[ct]]@result, file.path(FUN_DIR, paste0("GSEA_BP_", safe_ct, ".csv")))
}

cat("Wrote Mod7 figures/tables to", OUT_DIR, "\n")
