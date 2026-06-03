#!/usr/bin/env Rscript
# Talapas analysis pipeline 06 — parallels laptop Tutorial 06 (Pseudobulk DE, DESeq2).
# Learning notebook: Exercise_Folder/Tutorial_06_DESeq2_DE.qmd
# Run:  sbatch --job-name=pb_de --time=04:00:00 --mem=96G run_rscript.sbatch 06_pseudobulk_de.R
# In:   ../objects/nsclc_integrated.rds   Out: ../objects/nsclc_pseudobulk_de.csv
#
# NSCLC has no real biological condition, so we split the pseudo-donors into a
# synthetic two-group contrast purely to exercise the workflow. EXPECT ~no
# significant hits: this validates the plumbing, not biology.

suppressPackageStartupMessages({ library(Seurat); library(DESeq2); library(tidyverse) })
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
seu <- readRDS(file.path(OBJ_DIR, "nsclc_integrated.rds"))
DefaultAssay(seu) <- "RNA"

if (!"condition" %in% colnames(seu@meta.data))
  seu$condition <- ifelse(seu$donor %in% c("d1","d2","d3"), "GroupA", "GroupB")
seu$celltype <- seu$celltype_manual

# Step 3 — aggregate to (donor x condition x celltype) pseudobulk
pb <- AggregateExpression(seu, assays = "RNA", layer = "counts",
                          group.by = c("donor","condition","celltype"),
                          return.seurat = FALSE)$RNA

fix_id <- function(x) gsub("[ /]", "-", x)
group_meta <- seu@meta.data |>
  distinct(donor, condition, celltype) |>
  mutate(group_id = fix_id(paste(donor, condition, celltype, sep = "_")))
colnames(pb) <- fix_id(colnames(pb))
meta_pb <- tibble(group_id = colnames(pb)) |> left_join(group_meta, by = "group_id")

# Drop low-cell pseudobulk columns (< 10 cells)
cells_per_group <- seu@meta.data |>
  count(donor, condition, celltype) |>
  mutate(group_id = fix_id(paste(donor, condition, celltype, sep = "_")))
meta_pb <- meta_pb |>
  left_join(cells_per_group |> select(group_id, n_cells = n), by = "group_id")
keep <- !is.na(meta_pb$n_cells) & meta_pb$n_cells >= 10
pb <- pb[, keep]; meta_pb <- meta_pb[keep, ]

# Step 5 — DESeq2 per cell type
run_de <- function(ct) {
  k <- meta_pb$celltype == ct
  if (sum(k) < 4 || length(unique(meta_pb$condition[k])) < 2) return(NULL)
  dds <- DESeqDataSetFromMatrix(countData = pb[, k],
                                colData = as.data.frame(meta_pb[k, ]),
                                design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition","GroupB","GroupA"))
  res <- lfcShrink(dds, coef = "condition_GroupB_vs_GroupA", res = res, type = "apeglm")
  as.data.frame(res) |> rownames_to_column("gene") |> mutate(celltype = ct)
}
de_all <- map_dfr(unique(meta_pb$celltype), run_de) |> filter(!is.na(padj))

write_csv(de_all, file.path(OBJ_DIR, "nsclc_pseudobulk_de.csv"))
cat("Wrote", file.path(OBJ_DIR, "nsclc_pseudobulk_de.csv"), "with", nrow(de_all),
    "rows across", length(unique(de_all$celltype)), "cell types\n")
