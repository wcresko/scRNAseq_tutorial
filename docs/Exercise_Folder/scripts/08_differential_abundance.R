#!/usr/bin/env Rscript
# Talapas analysis pipeline 08 — parallels laptop Tutorial 08 (Differential Abundance, miloR).
# Learning notebook: Exercise_Folder/Tutorial_08_DifferentialAbundance.qmd
# Run:  sbatch --job-name=da --time=04:00:00 --mem=96G run_rscript.sbatch 08_differential_abundance.R
# In:   ../data/ifnb_integrated.rds   Out: ../data/ifnb_milo_da.csv, ../output/Mod8/
# Figures/tables (match the Mod8 notebook filenames): ../output/Mod8/Mod8_C*_*
#
# ifnb: real STIM vs CTRL contrast, 8 donors → genuine differential-abundance testing.

suppressPackageStartupMessages({
  library(Seurat); library(miloR); library(SingleCellExperiment)
  library(scater); library(tidyverse); library(patchwork)
})
set.seed(2026)
DATA_DIR <- Sys.getenv("DATA_DIR", "../data")
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod8") # figures/tables, named to match Tutorial_08.qmd
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("[dirs] data -> ", normalizePath(DATA_DIR), "  |  figures/tables -> ", normalizePath(OUT_DIR))

# Re-bind dplyr verbs masked by Bioconductor packages (S4Vectors, matrixStats, etc.).
select <- dplyr::select; filter <- dplyr::filter; count <- dplyr::count
desc   <- dplyr::desc;   slice  <- dplyr::slice
rename <- dplyr::rename; first  <- dplyr::first

# --- Figure saving --------------------------------------------------------
# save_fig(): write each figure as .png + .svg (vector). See Tutorial_08_DifferentialAbundance.qmd.
save_fig <- function(filename, plot, width, height, dpi = 300, ...) {
  ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, ...)
  svg_path <- paste0(tools::file_path_sans_ext(filename), ".svg")
  tryCatch(ggplot2::ggsave(svg_path, plot, width = width, height = height, ...),
           error = function(e) message("  (could not write ", basename(svg_path),
                                       " - install 'svglite'? ", conditionMessage(e), ")"))
}

seu <- readRDS(file.path(DATA_DIR, "ifnb_integrated.rds"))

# Real condition + replicate. miloR needs one sample ID per replicate: donor x stim.
if (!"donor" %in% colnames(seu@meta.data)) seu$donor <- seu$ind
seu$condition <- seu$stim
seu$sample    <- paste(seu$donor, seu$stim, sep = "_")

# Table out: cells per sample x condition (Mod8_C2)
as.data.frame(table(sample = seu$sample, stim = seu$stim)) |>
  rename(n_cells = Freq) |>
  write_csv(file.path(OUT_DIR, "Mod8_C2_cells_per_sample.csv"))

# Step 2 — convert to Milo
sce  <- as.SingleCellExperiment(seu)
milo <- Milo(sce)

# Step 3 — build kNN graph + neighbourhoods on the harmony embedding.
# k = 30, d = 20 are reasonable defaults for ~14k cells (kept in sync with
# the laptop notebook, Tutorial_08).
milo <- buildGraph(milo, k = 30, d = 20, reduced.dim = "HARMONY")
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 20, refined = TRUE,
                   reduced_dims = "HARMONY")

# Figure out: neighbourhood-size distribution (Mod8_C4)
p_nhood_size <- plotNhoodSizeHist(milo) +
  labs(title    = "Neighbourhood size distribution",
       subtitle = "Most neighbourhoods should hold 50-200 cells for a powered NB test",
       x        = "Neighbourhood size (cells)",
       y        = "Number of neighbourhoods",
       caption  = "Module 8 · Differential Abundance")
save_fig(file.path(OUT_DIR, "Mod8_C4_nhood_size_hist.png"), p_nhood_size,
       width = 6, height = 4, dpi = 300)

# Step 4 — count cells per neighbourhood per sample
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample = "sample")

# Step 5 — design + test
design_df <- as.data.frame(colData(milo)) |>
  distinct(sample, condition) |> column_to_rownames("sample")
milo <- calcNhoodDistance(milo, d = 20, reduced.dim = "HARMONY")
res  <- testNhoods(milo, design = ~ condition, design.df = design_df)

# Table out: top-10 neighbourhoods by SpatialFDR (Mod8_C6)
res |>
  arrange(SpatialFDR) |>
  head(10) |>
  write_csv(file.path(OUT_DIR, "Mod8_C6_da_top10.csv"))

# Step 6 — annotate neighbourhoods with the majority cell-type label
res <- annotateNhoods(milo, res, coldata_col = "seurat_annotations")

# Table out: DA neighbourhood counts per cell type, by direction (Mod8_C7)
res |>
  filter(SpatialFDR < 0.1) |>
  count(seurat_annotations, sign = sign(logFC)) |>
  pivot_wider(names_from = sign, values_from = n, names_prefix = "logFC_") |>
  rename(up_in_STIM = logFC_1, up_in_CTRL = `logFC_-1`) |>
  write_csv(file.path(OUT_DIR, "Mod8_C7_da_by_celltype.csv"))

# Step 7 — neighbourhood graph: cell-type UMAP alongside the DA graph (Mod8_C8)
milo <- buildNhoodGraph(milo)
um <- plotReducedDim(milo, dimred = "UMAP_HARMONY",
                     colour_by = "seurat_annotations",
                     text_by   = "seurat_annotations") +
  ggtitle("Cell types") +
  guides(colour = "none")
da <- plotNhoodGraphDA(milo, res, alpha = 0.1, layout = "UMAP_HARMONY") +
  ggtitle("DA neighbourhoods (red = up STIM, blue = up CTRL)")
p_graph <- (um + da) +
  plot_annotation(
    title    = "Differential abundance across the harmony UMAP — STIM vs CTRL",
    subtitle = "Each node is a neighbourhood, sized by cell count, coloured by log-fold change",
    caption  = "Module 8 · Differential Abundance")
save_fig(file.path(OUT_DIR, "Mod8_C8_da_nhood_graph.png"), p_graph,
       width = 12, height = 5, dpi = 300)

# Step 8 — beeswarm: per-neighbourhood logFC grouped by cell type (Mod8_C9)
p_beeswarm <- plotDAbeeswarm(res, group.by = "seurat_annotations") +
  labs(title    = "Differential abundance per cell type — STIM vs CTRL",
       subtitle = "Each dot is a neighbourhood; above zero = enriched in STIM, below = enriched in CTRL",
       x        = "Cell type",
       y        = "Log-fold change (STIM vs CTRL)",
       colour   = "SpatialFDR")
save_fig(file.path(OUT_DIR, "Mod8_C9_da_beeswarm.png"), p_beeswarm,
       width = 8, height = 6, dpi = 300)

# Pipeline hand-off: keep the DATA_DIR copy so downstream steps pick it up unchanged.
write_csv(as_tibble(res), file.path(DATA_DIR, "ifnb_milo_da.csv"))

# Table out: full DA result table for this module's outputs (Mod8_C10)
write_csv(as_tibble(res), file.path(OUT_DIR, "Mod8_C10_milo_da_results.csv"))

cat("Wrote", file.path(DATA_DIR, "ifnb_milo_da.csv"), "\n")
cat("Wrote Mod8 figures/tables to", OUT_DIR, "\n")
