#!/usr/bin/env Rscript
# Talapas analysis pipeline 15 — parallels laptop Tutorial 15 (scATAC-seq with Signac).
# Standalone bonus track. Dataset: PBMC 10k scATAC-seq v1 (hg19) in ../data/.
# Learning notebook: Exercise_Folder/Tutorial_15_scATACseq_Signac.qmd
# Run:  sbatch --job-name=atac --mem=64G run_rscript.sbatch 15_scatac_signac.R
# Out:  ../data/pbmc_atac_clustered.rds, ../output/Mod15/
# Figures/tables (match the Mod15 notebook filenames): ../output/Mod15/Mod15_C*_*

suppressPackageStartupMessages({
  library(Signac); library(Seurat); library(EnsDb.Hsapiens.v75)
  library(GenomeInfoDb); library(tidyverse); library(patchwork)
})
set.seed(2026)
DATA_DIR <- Sys.getenv("DATA_DIR", "../data")
OUT_DIR  <- Sys.getenv("OUT_DIR",  "../output/Mod15") # figures/tables, named to match Tutorial_15.qmd
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)  # input datasets land here
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("[dirs] data -> ", normalizePath(DATA_DIR), "  |  figures/tables -> ", normalizePath(OUT_DIR))

# Step 1 — Peek at the raw fragment file
frag.file <- read.delim(file.path(DATA_DIR, "atac_v1_pbmc_10k_fragments.tsv.gz"),
                        header = FALSE, nrows = 10)

# Table out: first rows of the raw fragment file (Mod15_C2)
frag_preview <- frag.file
colnames(frag_preview) <- c("chrom", "start", "end", "cell_barcode", "duplicate_count")
write_csv(as_tibble(frag_preview), file.path(OUT_DIR, "Mod15_C2_fragment_preview.csv"))

# Step 2 — Load peaks matrix + create the ChromatinAssay object
counts <- Read10X_h5(file.path(DATA_DIR, "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"))
chrom_assay <- CreateChromatinAssay(
  counts = counts, sep = c(":", "-"),
  fragments = file.path(DATA_DIR, "atac_v1_pbmc_10k_fragments.tsv.gz"),
  min.cells = 10, min.features = 200)
metadata <- read.csv(file.path(DATA_DIR, "atac_v1_pbmc_10k_singlecell.csv"),
                     header = TRUE, row.names = 1)
pbmc <- CreateSeuratObject(counts = chrom_assay, meta.data = metadata, assay = "ATAC")

# Step 3 — Gene annotation (hg19 -> EnsDb v75; convert to UCSC chr naming)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
Annotation(pbmc) <- annotations

# Step 4 — ATAC-specific QC metrics
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc, fast = FALSE)
pbmc$blacklist_ratio    <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

# Step 5 — Visualize QC (on the un-filtered object, matching the notebook) (Mod15_C6)
a1 <- DensityScatter(pbmc, x = "nCount_ATAC", y = "TSS.enrichment",
                     log_x = TRUE, quantiles = TRUE) +
  labs(title = "Depth vs TSS enrichment",
       x     = "nCount_ATAC (fragments per cell, log scale)",
       y     = "TSS enrichment score")
a2 <- DensityScatter(pbmc, x = "nucleosome_signal", y = "TSS.enrichment",
                     log_x = TRUE, quantiles = TRUE) +
  labs(title = "Nucleosome signal vs TSS enrichment",
       x     = "Nucleosome signal (log scale)",
       y     = "TSS enrichment score")
p_qc_density <- (a1 | a2) +
  plot_annotation(
    title    = "scATAC-seq QC density scatters — PBMC 10k v1",
    subtitle = "Good cells: sufficient depth, high TSS enrichment, low nucleosome signal",
    caption  = "Module 15 · scATAC-seq with Signac")
ggsave(file.path(OUT_DIR, "Mod15_C6_qc_density_scatter.png"), p_qc_density,
       width = 10, height = 5, dpi = 300)

p_qc_vln <- VlnPlot(pbmc,
        features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment",
                     "nucleosome_signal", "blacklist_ratio", "pct_reads_in_peaks"),
        pt.size = 0.1, ncol = 6) &
  theme(plot.title = element_text(size = 10))
p_qc_vln <- p_qc_vln +
  plot_annotation(
    title    = "Per-cell ATAC QC metric distributions — PBMC 10k v1",
    subtitle = "Depth, complexity, TSS enrichment, nucleosome signal, blacklist ratio and peak-read fraction",
    caption  = "Module 15 · scATAC-seq with Signac")
ggsave(file.path(OUT_DIR, "Mod15_C6_qc_violins.png"), p_qc_vln,
       width = 14, height = 4, dpi = 300)

# Step 6 — Filter (Signac PBMC vignette thresholds; tune per dataset)
pbmc <- subset(pbmc, subset = nCount_ATAC > 3000 & nCount_ATAC < 30000 &
                              pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 &
                              nucleosome_signal < 4 & TSS.enrichment > 3)

# Step 7 — TF-IDF normalization, top features, SVD (LSI)
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)

# Figure out: LSI component vs sequencing-depth correlation (Mod15_C8)
p_depthcor <- DepthCor(pbmc) +
  labs(title    = "LSI component correlation with sequencing depth",
       subtitle = "Component 1 typically tracks depth, not biology — exclude it downstream",
       x        = "LSI component",
       y        = "Correlation with depth")
ggsave(file.path(OUT_DIR, "Mod15_C8_depthcor.png"), p_depthcor,
       width = 7, height = 4, dpi = 300)

# Step 8 — UMAP + clustering on LSI dims 2:30 (drop depth-correlated component 1)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindClusters(pbmc, algorithm = 3)   # SLM

# Figure out: UMAP coloured by SLM cluster (Mod15_C9)
p_umap <- DimPlot(pbmc, label = TRUE) +
  NoLegend() +
  labs(title    = "scATAC-seq clusters — PBMC 10k v1",
       subtitle = "UMAP on LSI components 2:30, clustered with SLM (algorithm 3)",
       x        = "UMAP 1",
       y        = "UMAP 2")
ggsave(file.path(OUT_DIR, "Mod15_C9_umap_clusters.png"), p_umap,
       width = 7, height = 6, dpi = 300)

saveRDS(pbmc, file.path(DATA_DIR, "pbmc_atac_clustered.rds"))
cat("Wrote", file.path(DATA_DIR, "pbmc_atac_clustered.rds"), "with", ncol(pbmc), "cells\n")
cat("Wrote Mod15 figures/tables to", OUT_DIR, "\n")
