#!/usr/bin/env Rscript
# Talapas 07 — Pseudobulk DE with DESeq2
# Aggregates single-cell counts to (donor x condition x cell type), then runs
# DESeq2 per cell type. Writes nsclc_pseudobulk_de.csv.
#
# The NSCLC dataset has 7 donors but no biological condition variable, so this
# script demonstrates the workflow against a synthetic two-condition split.
# Replace the `condition` assignment with your real metadata for real DE.

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(tidyverse)
})

set.seed(2026)

OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "nsclc_integrated.rds"))
DefaultAssay(seu) <- "RNA"

# Synthetic condition for teaching: split donors d1-d3 vs d4-d7
if (!"condition" %in% colnames(seu@meta.data)) {
  seu$condition <- ifelse(seu$donor %in% c("d1","d2","d3"), "GroupA", "GroupB")
}
seu$celltype <- seu$celltype_manual

# Aggregate
pb <- AggregateExpression(
  seu, assays = "RNA", slot = "counts",
  group.by = c("donor", "condition", "celltype"),
  return.seurat = FALSE)$RNA

# Build metadata table for the pseudobulk columns
fix_id <- function(x) gsub("[ /]", "-", x)
group_meta <- seu@meta.data |>
  distinct(donor, condition, celltype) |>
  mutate(group_id = fix_id(paste(donor, condition, celltype, sep = "_")))
colnames(pb) <- fix_id(colnames(pb))
meta_pb <- tibble(group_id = colnames(pb)) |>
  left_join(group_meta, by = "group_id")

# Filter low-cell groups
cells_per_group <- seu@meta.data |>
  count(donor, condition, celltype) |>
  mutate(group_id = fix_id(paste(donor, condition, celltype, sep = "_")))
meta_pb <- meta_pb |>
  left_join(cells_per_group |> select(group_id, n_cells = n), by = "group_id")
keep <- !is.na(meta_pb$n_cells) & meta_pb$n_cells >= 10
pb <- pb[, keep]
meta_pb <- meta_pb[keep, ]

run_de <- function(ct) {
  k <- meta_pb$celltype == ct
  if (sum(k) < 4 || length(unique(meta_pb$condition[k])) < 2) return(NULL)
  dds <- DESeqDataSetFromMatrix(
    countData = pb[, k],
    colData   = as.data.frame(meta_pb[k, ]),
    design    = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "GroupB", "GroupA"))
  res <- lfcShrink(dds, coef = "condition_GroupB_vs_GroupA",
                   res = res, type = "apeglm")
  as.data.frame(res) |>
    rownames_to_column("gene") |>
    mutate(celltype = ct)
}

de_all <- map_dfr(unique(meta_pb$celltype), run_de) |> filter(!is.na(padj))
write_csv(de_all, file.path(OBJ_DIR, "nsclc_pseudobulk_de.csv"))
cat("Wrote", file.path(OBJ_DIR, "nsclc_pseudobulk_de.csv"),
    "with", nrow(de_all), "rows across",
    length(unique(de_all$celltype)), "cell types\n")
