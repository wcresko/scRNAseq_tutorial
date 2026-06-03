#!/usr/bin/env Rscript
# Talapas analysis pipeline 10 — parallels laptop Tutorial 10 (WGCNA).
# Standalone bonus track. Dataset: GSE152418 COVID-19 PBMC bulk RNA-seq.
# Learning notebook: Exercise_Folder/Tutorial_10_WGCNA.qmd
# Run:  sbatch --job-name=wgcna --mem=64G run_rscript.sbatch 10_wgcna.R
# In:   ../data/GSE152418_p20047_Study1_RawCounts.txt   Out: ../objects/wgcna_*.csv

suppressPackageStartupMessages({
  library(WGCNA); library(DESeq2); library(GEOquery); library(tidyverse)
})
set.seed(2026)
DATA_DIR <- Sys.getenv("DATA_DIR", "../data")
OBJ_DIR  <- Sys.getenv("OBJ_DIR",  "../objects")
dir.create(OBJ_DIR, showWarnings = FALSE, recursive = TRUE)

# Step 1 — Load counts + GEO metadata, reshape to a counts matrix
data <- read.delim(file.path(DATA_DIR, "GSE152418_p20047_Study1_RawCounts.txt"), header = TRUE)
gse  <- getGEO("GSE152418", GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))[, c(1, 2, 46:50)]
data <- data |>
  gather(key = "samples", value = "counts", -ENSEMBLID) |>
  mutate(samples = gsub("\\.", "-", samples)) |>
  inner_join(phenoData, by = c("samples" = "title")) |>
  select(1, 3, 4) |>
  spread(key = "geo_accession", value = "counts") |>
  column_to_rownames("ENSEMBLID")

# Step 2 — Outlier detection; drop bad genes + obvious outlier samples
gsg  <- goodSamplesGenes(t(data)); data <- data[gsg$goodGenes == TRUE, ]
samples.to.exclude <- c("GSM4614993", "GSM4614994", "GSM4614995")
data.subset <- data[, !(colnames(data) %in% samples.to.exclude)]
colData <- phenoData |> filter(!(geo_accession %in% samples.to.exclude)) |>
  column_to_rownames("geo_accession")

# Step 3 — DESeq2 VST (unsupervised network -> design ~ 1)
dds   <- DESeqDataSetFromMatrix(countData = data.subset, colData = colData, design = ~ 1)
dds75 <- dds[rowSums(counts(dds) >= 15) >= 24, ]
norm.counts <- assay(vst(dds75)) |> t()    # WGCNA wants samples x genes

# Step 4 — Pick soft-thresholding power (inspect the table, pick smallest with R^2 ~ 0.8)
sft <- pickSoftThreshold(norm.counts, powerVector = c(1:10, seq(12, 30, by = 2)),
                         networkType = "signed", verbose = 5)
write_csv(sft$fitIndices, file.path(OBJ_DIR, "wgcna_softpower.csv"))

# Step 5 — Build the network (beta = 18 from Patel's example for this dataset)
soft_power <- 18
temp_cor <- cor; cor <- WGCNA::cor          # WGCNA shadows base::cor
bwnet <- blockwiseModules(norm.counts, maxBlockSize = 14000, TOMType = "signed",
                          power = soft_power, mergeCutHeight = 0.25,
                          numericLabels = FALSE, randomSeed = 1234, verbose = 3)
cor <- temp_cor                              # restore base::cor

# Step 6 — Module eigengenes + module-trait correlation
module_eigengenes <- bwnet$MEs
traits <- binarizeCategoricalColumns(colData$`severity:ch1`, includePairwise = FALSE,
                                     includeLevelVsAll = TRUE, minCount = 1)
mtc <- cor(module_eigengenes, traits, use = "p")
write_csv(as.data.frame(mtc) |> rownames_to_column("module"),
          file.path(OBJ_DIR, "wgcna_module_trait_cor.csv"))
write_csv(as.data.frame(bwnet$colors) |> rownames_to_column("gene") |>
            rename(module = `bwnet$colors`),
          file.path(OBJ_DIR, "wgcna_gene_modules.csv"))

cat("Wrote WGCNA soft-power, module assignments, and trait-correlation tables to", OBJ_DIR, "\n")
