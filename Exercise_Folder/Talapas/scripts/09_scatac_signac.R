#!/usr/bin/env Rscript
# Talapas analysis pipeline 09 — parallels laptop Tutorial 09 (scATAC-seq with Signac).
# Standalone bonus track. Dataset: PBMC 10k scATAC-seq v1 (hg19) in ../data/.
# Learning notebook: Exercise_Folder/Tutorial_09_scATACseq_Signac.qmd
# Run:  sbatch --job-name=atac --mem=64G run_rscript.sbatch 09_scatac_signac.R
# Out:  ../objects/pbmc_atac_clustered.rds

suppressPackageStartupMessages({
  library(Signac); library(Seurat); library(EnsDb.Hsapiens.v75)
  library(GenomeInfoDb); library(tidyverse)
})
set.seed(2026)
DATA_DIR <- Sys.getenv("DATA_DIR", "../data")
OBJ_DIR  <- Sys.getenv("OBJ_DIR",  "../objects")
dir.create(OBJ_DIR, showWarnings = FALSE, recursive = TRUE)

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

# Step 6 — Filter (Signac PBMC vignette thresholds; tune per dataset)
pbmc <- subset(pbmc, subset = nCount_ATAC > 3000 & nCount_ATAC < 30000 &
                              pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 &
                              nucleosome_signal < 4 & TSS.enrichment > 3)

# Step 7 — TF-IDF normalization, top features, SVD (LSI)
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)

# Step 8 — UMAP + clustering on LSI dims 2:30 (drop depth-correlated component 1)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindClusters(pbmc, algorithm = 3)   # SLM

saveRDS(pbmc, file.path(OBJ_DIR, "pbmc_atac_clustered.rds"))
cat("Wrote", file.path(OBJ_DIR, "pbmc_atac_clustered.rds"), "with", ncol(pbmc), "cells\n")
