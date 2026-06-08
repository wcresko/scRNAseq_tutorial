#!/usr/bin/env Rscript
# Talapas analysis pipeline 06 — parallels laptop Tutorial 06 (Pseudobulk DE, DESeq2).
# Learning notebook: Exercise_Folder/Tutorial_06_DESeq2_DE.qmd
# Run:  sbatch --job-name=pb_de --time=04:00:00 --mem=96G run_rscript.sbatch 06_pseudobulk_de.R
# In:   ../data/ifnb_integrated.rds
# Out:  ../data/ifnb_pseudobulk_de.csv  +  ../output/Mod6 figures/tables (match Tutorial_06.qmd)
#
# ifnb has a REAL biological contrast (IFN-beta STIM vs CTRL) and REAL donors
# (8 lupus patients, `ind`), so this is genuine pseudobulk DE — expect a strong
# interferon-stimulated-gene (ISG) signature, not a null result.

suppressPackageStartupMessages({ library(Seurat); library(DESeq2); library(tidyverse); library(patchwork) })
set.seed(2026)
DATA_DIR <- Sys.getenv("DATA_DIR", "../data")
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod6") # figures/tables, named to match Tutorial_06.qmd
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)   # pipeline hand-off objects (.rds/.csv)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
message("[dirs] data -> ", normalizePath(DATA_DIR), "  |  figures/tables -> ", normalizePath(OUT_DIR))
seu <- readRDS(file.path(DATA_DIR, "ifnb_integrated.rds"))
DefaultAssay(seu) <- "RNA"

# Real donor (ind) and real condition (stim)
if (!"donor" %in% colnames(seu@meta.data)) seu$donor <- seu$ind
seu$condition <- seu$stim                 # CTRL / STIM
seu$celltype  <- seu$seurat_annotations

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

# Step 5 — DESeq2 per cell type (STIM vs CTRL)
run_de <- function(ct) {
  k <- meta_pb$celltype == ct
  if (sum(k) < 4 || length(unique(meta_pb$condition[k])) < 2) return(NULL)
  cd <- as.data.frame(meta_pb[k, ])
  cd$condition <- factor(cd$condition, levels = c("CTRL", "STIM"))
  dds <- DESeqDataSetFromMatrix(countData = pb[, k], colData = cd, design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- lfcShrink(dds, coef = "condition_STIM_vs_CTRL", type = "apeglm")
  as.data.frame(res) |> rownames_to_column("gene") |> mutate(celltype = ct)
}
de_all <- map_dfr(unique(meta_pb$celltype), run_de) |> filter(!is.na(padj))

# Pipeline hand-off: script 07 reads this from DATA_DIR — keep it there.
write_csv(de_all, file.path(DATA_DIR, "ifnb_pseudobulk_de.csv"))
cat("Wrote", file.path(DATA_DIR, "ifnb_pseudobulk_de.csv"), "with", nrow(de_all),
    "rows across", length(unique(de_all$celltype)), "cell types\n")

# ===========================================================================
# Mod6 figures/tables — reuse the Tutorial_06.qmd filenames (ifnb pseudobulk
# section, chunks C16-C21). The qmd's airway teaching section (C3-C11) is not
# part of this pipeline and is intentionally not reproduced here.
# ===========================================================================

# --- Table: pseudobulk column metadata, surviving groups (Mod6_C16) --------
write_csv(meta_pb, file.path(OUT_DIR, "Mod6_C16_pseudobulk_coldata.csv"))

# --- Table: surviving pseudobulk groups per celltype x condition (Mod6_C17) -
as.data.frame(table(celltype = meta_pb$celltype, condition = meta_pb$condition)) |>
  write_csv(file.path(OUT_DIR, "Mod6_C17_groups_per_celltype.csv"))

# --- Table: full per-cell-type pseudobulk DE results (Mod6_C18) ------------
write_csv(de_all, file.path(OUT_DIR, "Mod6_C18_pseudobulk_de_all.csv"))

# --- Table: top 10 hits per cell type by padj (Mod6_C19) -------------------
de_top <- de_all |>
  group_by(celltype) |>
  slice_min(padj, n = 10) |>
  ungroup()
write_csv(de_top, file.path(OUT_DIR, "Mod6_C19_pseudobulk_de_top.csv"))

# --- Figure: per-cell-type pseudobulk volcano plots (Mod6_C20) -------------
plot_types <- intersect(c("CD14 Mono", "CD4 Naive T", "B"), unique(de_all$celltype))
p_volcano <- de_all |>
  filter(celltype %in% plot_types) |>
  ggplot(aes(log2FoldChange, -log10(padj))) +
  geom_point(alpha = 0.4, size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = 3, color = "grey50") +
  facet_wrap(~ celltype) +
  labs(title    = "Pseudobulk DE: STIM vs CTRL (log2FC > 0 = up in STIM)",
       subtitle = "Dashed line: padj = 0.05; dotted lines: |log2FC| = 1",
       x        = "log2 fold change (apeglm-shrunken)",
       y        = "-log10(padj)",
       caption  = "Module 6 · DESeq2 pseudobulk DE")
ggsave(file.path(OUT_DIR, "Mod6_C20_pseudobulk_volcano.png"), p_volcano,
       width = 10, height = 4, dpi = 300)

# --- Tables: student-facing DE copies (Mod6_C21) ---------------------------
write_csv(de_all, file.path(OUT_DIR, "Mod6_C21_pseudobulk_de.csv"))
write_csv(de_top, file.path(OUT_DIR, "Mod6_C21_pseudobulk_de_top.csv"))

cat("Wrote Mod6 figures/tables to", OUT_DIR, "\n")
