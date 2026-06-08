#!/usr/bin/env Rscript
# Talapas analysis pipeline 01 — parallels laptop Tutorial 01 (QC & Preprocessing).
# Dataset: ifnb (Kang et al. 2017) — the SAME dataset as the laptop tutorials.
# Learning notebook: Exercise_Folder/Tutorial_01_QC_Preprocessing.qmd  (Module 1)
# Run:  sbatch --job-name=qc run_rscript.sbatch 01_qc_preprocessing.R
# In:   muscData::Kang18_8vs8() (ExperimentHub cache)   Out: ../data/ifnb_preprocessed.rds
# Figures/tables (match the Mod1 notebook filenames): ../output/Mod1/Mod1_C*_*

suppressPackageStartupMessages({
  library(Seurat); library(muscData); library(SingleCellExperiment)
  library(tidyverse); library(patchwork)
})
set.seed(2026)

DATA_DIR <- Sys.getenv("DATA_DIR", "../data")   # every step reads/writes pipeline objects here
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod1") # figures/tables, named to match Tutorial_01.qmd
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
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

# Step 1 — Load ifnb from Bioconductor's ExperimentHub (cached after the first
# download). Keep singlets that have an author-assigned cell-type label, and mirror
# the metadata names the rest of the series expects:
#   - `stim`               CTRL / STIM condition
#   - `seurat_annotations` author-curated cell types
#   - `donor` (= `ind`)    the 8 real lupus donors used downstream for DE / DA
sce <- Kang18_8vs8()
sce <- sce[, sce$multiplets == "singlet" & !is.na(sce$cell)]
# muscData's Kang18_8vs8 stores features as SYMBOL_ENSEMBLID (e.g.
# "ISG15_ENSG00000187608"). Keep only the SYMBOL so markers / Azimuth /
# org.Hs.eg.db match. Strip ONLY a trailing "_ENSEMBLID" suffix (NOT everything
# after the first dash — real symbols like HLA-A / MT-CO1 contain dashes); keep
# the original name where there is no symbol; de-duplicate collapsed symbols.
cm   <- counts(sce)
orig <- rownames(cm)
sym  <- sub("_ENS[A-Z]*[0-9]+(\\.[0-9]+)?$", "", orig)
blank <- is.na(sym) | sym == "" | grepl("^ENS[A-Z]*[0-9]+", sym)
sym[blank] <- orig[blank]
rownames(cm) <- make.unique(gsub("_", "-", sym))
if (sum(rownames(cm) != orig) > 0)
  message(sum(rownames(cm) != orig), " feature names simplified to gene symbols")
seu <- CreateSeuratObject(counts    = cm,
                          meta.data = as.data.frame(colData(sce)),
                          project   = "ifnb", min.cells = 3, min.features = 200)
seu$stim               <- factor(toupper(seu$stim), levels = c("CTRL", "STIM"))
seu$seurat_annotations <- factor(seu$cell)
seu$donor              <- seu$ind   # real biological replicate

# Table out: cells per sample (Mod1_C2)
enframe(table(seu$stim), name = "sample", value = "n_cells") |>
  mutate(n_cells = as.integer(n_cells)) |>
  write_csv(file.path(OUT_DIR, "Mod1_C2_cells_per_sample.csv"))

# Step 3 — QC metrics (SoupX needs a raw matrix and is illustrative-only on ifnb;
# scDblFinder is optional — see the notebook).
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

# Figure out: per-sample QC violins (Mod1_C5). Plotted pre-normalization on QC
# metadata, so wrap to silence the expected "data layer empty, using counts" note.
p_qc_vln <- suppressWarnings(VlnPlot(seu,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
        group.by = "stim", ncol = 4, pt.size = 0) &
  xlab("Sample (stim condition)") &
  theme(plot.title = element_text(size = 11)))
p_qc_vln <- p_qc_vln + plot_annotation(
  title    = "Per-cell QC metrics by sample — ifnb (CTRL vs IFN-beta STIM)",
  subtitle = "Genes/cell, UMIs/cell, mitochondrial % and ribosomal % distributions",
  caption  = "Module 1 · QC & Preprocessing (Talapas pipeline)")
save_fig(file.path(OUT_DIR, "Mod1_C5_qc_violins.png"), p_qc_vln, width = 12, height = 4, dpi = 300)

# Table out: per-sample QC summary (Mod1_C4)
seu@meta.data |>
  group_by(stim) |>
  summarise(n_cells = n(),
            median_nFeature = median(nFeature_RNA),
            median_nCount   = median(nCount_RNA),
            median_percent_mt = round(median(percent.mt), 2),
            median_percent_rb = round(median(percent.rb), 2),
            .groups = "drop") |>
  write_csv(file.path(OUT_DIR, "Mod1_C4_qc_metric_summary.csv"))

# Step 7 — Filter (PBMC-appropriate thresholds)
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Table out: cells surviving the filter, per sample (Mod1_C8)
enframe(table(seu$stim), name = "sample", value = "n_cells_after_filter") |>
  mutate(n_cells_after_filter = as.integer(n_cells_after_filter)) |>
  write_csv(file.path(OUT_DIR, "Mod1_C8_cells_after_filter.csv"))

# Step 8 — Normalize
seu <- NormalizeData(seu)

# Step 9 — Highly variable features
seu <- FindVariableFeatures(seu, nfeatures = 2000)

# Figure + table out: HVGs (Mod1_C10)
top10 <- head(VariableFeatures(seu), 10)
p_hvg <- LabelPoints(plot = VariableFeaturePlot(seu), points = top10, repel = TRUE) +
  labs(title = "Highly variable genes (vst) — top 2,000 selected",
       subtitle = "Top 10 most variable genes labelled (many are interferon-stimulated genes)",
       x = "Average expression (log scale)", y = "Standardized variance",
       colour = "Selected as HVG")
save_fig(file.path(OUT_DIR, "Mod1_C10_variable_features.png"), p_hvg, width = 7, height = 5, dpi = 300)
tibble(rank = seq_along(top10), gene = top10) |>
  write_csv(file.path(OUT_DIR, "Mod1_C10_top10_hvgs.csv"))

# Step 10 — Scale
seu <- ScaleData(seu, vars.to.regress = c("nCount_RNA", "percent.mt"))

saveRDS(seu, file.path(DATA_DIR, "ifnb_preprocessed.rds"))
cat("Wrote", file.path(DATA_DIR, "ifnb_preprocessed.rds"), "with", ncol(seu), "cells\n")
cat("Wrote Mod1 figures/tables to", OUT_DIR, "\n")
