#!/usr/bin/env Rscript
# Talapas analysis pipeline 13 — parallels laptop Tutorial 13 (WGCNA).
# Standalone bonus track. Dataset: GSE152418 COVID-19 PBMC bulk RNA-seq.
# Learning notebook: Exercise_Folder/Tutorial_13_WGCNA.qmd
# Run:  sbatch --job-name=wgcna --mem=64G run_rscript.sbatch 13_wgcna.R
# In:   ../data/GSE152418_p20047_Study1_RawCounts.txt
# Out:  ../objects/wgcna_*.csv, ../output/Mod13/
# Figures/tables (match the Mod13 notebook filenames): ../output/Mod13/Mod13_C*_*

suppressPackageStartupMessages({
  library(WGCNA); library(DESeq2); library(GEOquery); library(tidyverse)
  library(CorLevelPlot); library(gridExtra); library(patchwork)
})
set.seed(2026)
DATA_DIR <- Sys.getenv("DATA_DIR", "../data")
OBJ_DIR  <- Sys.getenv("OBJ_DIR",  "../objects")
OUT_DIR  <- Sys.getenv("OUT_DIR",  "../output/Mod13") # figures/tables, named to match Tutorial_13.qmd
dir.create(OBJ_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Step 1 — Load counts + GEO metadata, reshape to a counts matrix
data <- read.delim(file.path(DATA_DIR, "GSE152418_p20047_Study1_RawCounts.txt"), header = TRUE)
gse  <- getGEO("GSE152418", GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))[, c(1, 2, 46:50)]

# Table out: the trimmed sample metadata (phenoData) (Mod13_C3)
phenoData |>
  rownames_to_column("geo_accession_rowname") |>
  write_csv(file.path(OUT_DIR, "Mod13_C3_pheno_data.csv"))

data <- data |>
  gather(key = "samples", value = "counts", -ENSEMBLID) |>
  mutate(samples = gsub("\\.", "-", samples)) |>
  inner_join(phenoData, by = c("samples" = "title")) |>
  select(1, 3, 4) |>
  spread(key = "geo_accession", value = "counts") |>
  column_to_rownames("ENSEMBLID")

# Table out: counts-matrix dimensions (genes x samples) (Mod13_C4)
tibble(stage   = "Step 1 — tidy counts matrix (Ensembl x GEO accession)",
       genes   = nrow(data),
       samples = ncol(data)) |>
  write_csv(file.path(OUT_DIR, "Mod13_C4_counts_dimensions.csv"))

# Step 2 — Outlier detection; drop bad genes + obvious outlier samples
gsg  <- goodSamplesGenes(t(data))

# Table out: goodSamplesGenes outcome (genes/samples flagged) (Mod13_C5)
tibble(all_ok         = gsg$allOK,
       n_good_genes   = sum(gsg$goodGenes),
       n_bad_genes    = sum(!gsg$goodGenes),
       n_good_samples = sum(gsg$goodSamples),
       n_bad_samples  = sum(!gsg$goodSamples)) |>
  write_csv(file.path(OUT_DIR, "Mod13_C5_good_samples_genes.csv"))

data <- data[gsg$goodGenes == TRUE, ]

# Figure out (base graphics): sample dendrogram for outlier detection (Mod13_C6)
htree <- hclust(dist(t(data)), method = "average")
png(file.path(OUT_DIR, "Mod13_C6_sample_dendrogram.png"),
    width = 10, height = 6, units = "in", res = 300)
plot(htree,
     main = "Sample clustering for outlier detection — GSE152418 PBMCs",
     sub  = "Average linkage, Euclidean distance on samples",
     xlab = "Sample")
dev.off()

# Figure out: sample-level PCA (PC1 vs PC2) for outlier detection (Mod13_C7)
pca       <- prcomp(t(data))
pca.dat   <- pca$x
pca.var   <- pca$sdev^2
pca.var.p <- round(pca.var / sum(pca.var) * 100, digits = 2)
p_pca <- ggplot(as.data.frame(pca.dat), aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat), nudge_y = 5) +
  labs(title    = "Sample PCA for outlier detection — GSE152418 PBMCs",
       subtitle = "Each point is a sample; outliers sit far from the main cloud",
       x = paste0("PC1: ", pca.var.p[1], "%"),
       y = paste0("PC2: ", pca.var.p[2], "%"),
       caption  = "Module 13 · WGCNA") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "Mod13_C7_sample_pca.png"), p_pca, width = 7, height = 6, dpi = 300)

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
sft.data <- sft$fitIndices
# Pipeline hand-off: keep the OBJ_DIR copy.
write_csv(sft.data, file.path(OBJ_DIR, "wgcna_softpower.csv"))

# Figure out: soft-threshold scan — scale-free fit + mean connectivity (Mod13_C10)
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() + geom_text(nudge_y = 0.05) +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(title = "Scale-free topology fit",
       x = "Power", y = "Scale-free topology fit, signed R²") +
  theme_classic()
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() + geom_text(nudge_y = 0.1) +
  labs(title = "Mean connectivity",
       x = "Power", y = "Mean connectivity") +
  theme_classic()
p_sft <- (a1 / a2) +
  plot_annotation(
    title    = "Soft-thresholding power selection — GSE152418 PBMCs",
    subtitle = "Smallest power crossing signed R² ≈ 0.8 with reasonable connectivity",
    caption  = "Module 13 · WGCNA")
ggsave(file.path(OUT_DIR, "Mod13_C10_soft_threshold.png"), p_sft, width = 7, height = 8, dpi = 300)

# Step 5 — Build the network (beta = 18 from Patel's example for this dataset)
soft_power <- 18
temp_cor <- cor; cor <- WGCNA::cor          # WGCNA shadows base::cor
bwnet <- blockwiseModules(norm.counts, maxBlockSize = 14000, TOMType = "signed",
                          power = soft_power, mergeCutHeight = 0.25,
                          numericLabels = FALSE, randomSeed = 1234, verbose = 3)
cor <- temp_cor                              # restore base::cor

# Figure out (base graphics): gene dendrogram + module colours (Mod13_C12)
png(file.path(OUT_DIR, "Mod13_C12_module_dendrogram.png"),
    width = 10, height = 6, units = "in", res = 300)
plotDendroAndColors(bwnet$dendrograms[[1]],
                    cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide     = TRUE,
                    hang         = 0.03,
                    guideHang    = 0.05,
                    main         = "Gene clustering dendrogram & module colours")
dev.off()

# Step 6 — Module eigengenes + module-trait correlation
module_eigengenes <- bwnet$MEs

# Tables out: module eigengenes and genes-per-module (Mod13_C13)
module_eigengenes |>
  rownames_to_column("sample") |>
  write_csv(file.path(OUT_DIR, "Mod13_C13_module_eigengenes.csv"))
enframe(table(bwnet$colors), name = "module", value = "n_genes") |>
  mutate(n_genes = as.integer(n_genes)) |>
  arrange(desc(n_genes)) |>
  write_csv(file.path(OUT_DIR, "Mod13_C13_genes_per_module.csv"))

# Build the numeric trait matrix: a severe/not-severe binary plus per-level binarization.
traits <- colData |>
  mutate(disease_state_bin = ifelse(grepl("severe", `disease state:ch1`,
                                          ignore.case = TRUE), 1, 0)) |>
  select(disease_state_bin)
severity.out <- binarizeCategoricalColumns(colData$`severity:ch1`,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE, minCount = 1)
traits <- cbind(traits, severity.out)

nSamples <- nrow(norm.counts)
nGenes   <- ncol(norm.counts)
module.trait.corr   <- cor(module_eigengenes, traits, use = "p")
module.trait.corr.p <- corPvalueStudent(module.trait.corr, nSamples)

# Figure out (base graphics): module-trait correlation heatmap (Mod13_C15)
heatmap.data <- merge(module_eigengenes, traits, by = "row.names")
heatmap.data <- column_to_rownames(heatmap.data, var = "Row.names")
png(file.path(OUT_DIR, "Mod13_C15_module_trait_heatmap.png"),
    width = 8, height = 8, units = "in", res = 300)
CorLevelPlot(heatmap.data,
             x    = names(traits),
             y    = names(module_eigengenes),
             col  = c("blue1", "skyblue", "white", "pink", "red"),
             main = "Module-trait correlations — GSE152418 severity")
dev.off()

# Tables out: module-trait correlations + p-values (Mod13_C15)
as.data.frame(module.trait.corr) |>
  rownames_to_column("module") |>
  write_csv(file.path(OUT_DIR, "Mod13_C15_module_trait_corr.csv"))
as.data.frame(module.trait.corr.p) |>
  rownames_to_column("module") |>
  write_csv(file.path(OUT_DIR, "Mod13_C15_module_trait_pvalues.csv"))

# Pipeline hand-off: keep the OBJ_DIR copies so downstream steps pick them up unchanged.
write_csv(as.data.frame(module.trait.corr) |> rownames_to_column("module"),
          file.path(OBJ_DIR, "wgcna_module_trait_cor.csv"))
write_csv(as.data.frame(bwnet$colors) |> rownames_to_column("gene") |>
            rename(module = `bwnet$colors`),
          file.path(OBJ_DIR, "wgcna_gene_modules.csv"))

# Step 7 — Hub genes within a module (kME = gene-module eigengene correlation)
module.gene.mapping <- as.data.frame(bwnet$colors)
mm <- as.data.frame(cor(norm.counts, module_eigengenes, use = "p"))

# Tables out: gene-module assignments + top turquoise kME hub genes (Mod13_C16)
module.gene.mapping |>
  rownames_to_column("gene") |>
  write_csv(file.path(OUT_DIR, "Mod13_C16_gene_module_assignments.csv"))
mm |>
  rownames_to_column("gene") |>
  arrange(desc(MEturquoise)) |>
  head(15) |>
  write_csv(file.path(OUT_DIR, "Mod13_C16_turquoise_hub_genes.csv"))

cat("Wrote WGCNA soft-power, module assignments, and trait-correlation tables to", OBJ_DIR, "\n")
cat("Wrote Mod13 figures/tables to", OUT_DIR, "\n")
