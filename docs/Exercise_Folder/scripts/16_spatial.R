#!/usr/bin/env Rscript
# Talapas analysis pipeline 16 — parallels laptop Tutorial 16 (Spatial / Visium).
# Standalone bonus track. Dataset: stxBrain (10x Visium adult mouse brain) via SeuratData.
# Learning notebook: Exercise_Folder/Tutorial_16_Spatial_Transcriptomics.qmd
# Run:  sbatch --job-name=spatial --mem=64G run_rscript.sbatch 16_spatial.R
# Out:  ../objects/brain_spatial_integrated.rds
# Figures/tables (match the Mod16 notebook filenames): ../output/Mod16/Mod16_C*_*

suppressPackageStartupMessages({
  library(Seurat); library(SeuratData); library(tidyverse); library(patchwork)
})
set.seed(2026)
OBJ_DIR <- Sys.getenv("OBJ_DIR", "../objects")
OUT_DIR <- Sys.getenv("OUT_DIR", "../output/Mod16") # figures/tables, named to match Tutorial_16.qmd
dir.create(OBJ_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Step 1 — Load both sections (run SeuratData::InstallData("stxBrain") once first)
ant  <- LoadData("stxBrain", type = "anterior1")
post <- LoadData("stxBrain", type = "posterior1")

# Table out: per-section dimensions (features x spots) (Mod16_C2)
tibble(section  = c("anterior1", "posterior1"),
       features = c(nrow(ant), nrow(post)),
       spots    = c(ncol(ant), ncol(post))) |>
  write_csv(file.path(OUT_DIR, "Mod16_C2_section_dimensions.csv"))

# Step 2 — Spatial QC metric (mouse mito prefix "^mt-")
ant[["percent.mt"]]  <- PercentageFeatureSet(ant,  pattern = "^mt-")
post[["percent.mt"]] <- PercentageFeatureSet(post, pattern = "^mt-")

# Figure out: per-spot QC violins for the anterior section (Mod16_C3)
p_qc_vln <- VlnPlot(ant, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
                    pt.size = 0, ncol = 3) &
  xlab("Anterior section") &
  theme(plot.title = element_text(size = 11))
p_qc_vln <- p_qc_vln + plot_annotation(
  title    = "Per-spot QC metrics — stxBrain anterior section",
  subtitle = "UMIs/spot, genes/spot and mitochondrial % distributions",
  caption  = "Module 16 · Spatial Transcriptomics")
ggsave(file.path(OUT_DIR, "Mod16_C3_qc_violins.png"), p_qc_vln, width = 10, height = 4, dpi = 300)

# Figure out: spatial QC overlay on the tissue (Mod16_C3)
p_qc_spatial <- SpatialFeaturePlot(ant, features = c("nCount_Spatial", "percent.mt"),
                                   pt.size.factor = 1.6) + plot_layout(ncol = 2)
p_qc_spatial <- p_qc_spatial + plot_annotation(
  title    = "Spatial QC overlay — stxBrain anterior section",
  subtitle = "UMIs/spot and mitochondrial % on the histology image",
  caption  = "Module 16 · Spatial Transcriptomics")
ggsave(file.path(OUT_DIR, "Mod16_C3_qc_spatial.png"), p_qc_spatial, width = 10, height = 5, dpi = 300)

# Step 3 — SCTransform normalization (recommended for Visium's wider dynamic range)
ant  <- SCTransform(ant,  assay = "Spatial", verbose = FALSE)
post <- SCTransform(post, assay = "Spatial", verbose = FALSE)

# Step 4 — Figure out: known anatomical marker genes on the tissue (Mod16_C5)
p_markers_ant <- SpatialFeaturePlot(ant, features = c("Hpca", "Plp1", "Sst"),
                                    ncol = 3, alpha = c(0.1, 1)) +
  plot_annotation(
    title    = "Anatomical markers on tissue — anterior section",
    subtitle = "Hpca (hippocampus), Plp1 (white matter), Sst (interneurons)",
    caption  = "Module 16 · Spatial Transcriptomics")
ggsave(file.path(OUT_DIR, "Mod16_C5_markers_anterior.png"), p_markers_ant, width = 12, height = 4, dpi = 300)

p_markers_post <- SpatialFeaturePlot(post, features = c("Hpca", "Plp1", "Pcp4"),
                                     ncol = 3, alpha = c(0.1, 1)) +
  plot_annotation(
    title    = "Anatomical markers on tissue — posterior section",
    subtitle = "Hpca (hippocampus), Plp1 (white matter), Pcp4 (cerebellum)",
    caption  = "Module 16 · Spatial Transcriptomics")
ggsave(file.path(OUT_DIR, "Mod16_C5_markers_posterior.png"), p_markers_post, width = 12, height = 4, dpi = 300)

# Step 5 — Cluster spots on the anterior section (notebook Mod16_C6); also needed so the
# spatially variable features below run on a clustered object as in the notebook.
ant <- ant |>
  RunPCA(assay = "SCT", verbose = FALSE) |>
  FindNeighbors(reduction = "pca", dims = 1:30) |>
  FindClusters(resolution = 0.5, verbose = FALSE) |>
  RunUMAP(reduction = "pca", dims = 1:30)

# Figure out: anterior spot clusters in UMAP space and on tissue (Mod16_C6)
p_umap <- DimPlot(ant, reduction = "umap", label = TRUE) +
  labs(title = "Anterior — UMAP", x = "UMAP 1", y = "UMAP 2", colour = "Cluster")
p_tis  <- SpatialDimPlot(ant, label = TRUE, label.size = 3) +
  labs(title = "Anterior — clusters on tissue", fill = "Cluster")
p_cluster <- (p_umap + p_tis) + plot_annotation(
  title    = "Spot clustering — stxBrain anterior section",
  subtitle = "SNN clusters on the SCT assay, in UMAP space and on the histology image",
  caption  = "Module 16 · Spatial Transcriptomics")
ggsave(file.path(OUT_DIR, "Mod16_C6_anterior_clusters.png"), p_cluster, width = 12, height = 5, dpi = 300)

# Step 6 — Spatially variable features on the anterior section (Moran's I)
ant <- FindSpatiallyVariableFeatures(ant, assay = "SCT",
         features = VariableFeatures(ant)[1:1000], selection.method = "moransi")

# Figure + table out: top spatially variable genes (Mod16_C7)
top_sv <- SpatiallyVariableFeatures(ant, selection.method = "moransi")[1:6]
tibble(rank = seq_along(top_sv), gene = top_sv) |>
  write_csv(file.path(OUT_DIR, "Mod16_C7_spatially_variable_features.csv"))
p_svf <- SpatialFeaturePlot(ant, features = top_sv, ncol = 3, alpha = c(0.1, 1)) +
  plot_annotation(
    title    = "Top spatially variable features (Moran's I) — anterior section",
    subtitle = "Genes whose expression has the strongest spatial structure",
    caption  = "Module 16 · Spatial Transcriptomics")
ggsave(file.path(OUT_DIR, "Mod16_C7_spatially_variable_features.png"), p_svf, width = 12, height = 7, dpi = 300)

# Step 7 — Merge anterior + posterior and joint-cluster (SCT integration)
ant$slice <- "anterior"; post$slice <- "posterior"
brain.list <- list(ant = ant, post = post)
features   <- SelectIntegrationFeatures(brain.list, nfeatures = 3000, verbose = FALSE)
brain.list <- PrepSCTIntegration(brain.list, anchor.features = features, verbose = FALSE)
anchors    <- FindIntegrationAnchors(brain.list, normalization.method = "SCT",
                                     anchor.features = features, verbose = FALSE)
brain      <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
DefaultAssay(brain) <- "integrated"
brain <- brain |> RunPCA(verbose = FALSE) |> FindNeighbors(dims = 1:30) |>
  FindClusters(resolution = 0.5, verbose = FALSE) |> RunUMAP(dims = 1:30)

# Figure out: integrated UMAP coloured by section and joint cluster (Mod16_C8)
p_int_umap <- DimPlot(brain, group.by = c("slice", "seurat_clusters")) +
  plot_annotation(
    title    = "Integrated UMAP — anterior + posterior sections",
    subtitle = "Coloured by section of origin and by joint cluster",
    caption  = "Module 16 · Spatial Transcriptomics")
ggsave(file.path(OUT_DIR, "Mod16_C8_integrated_umap.png"), p_int_umap, width = 12, height = 5, dpi = 300)

# Figure out: joint clusters mapped back onto each section (Mod16_C8)
p_int_spatial <- SpatialDimPlot(brain, label = TRUE) +
  plot_annotation(
    title    = "Integrated clusters on tissue — both sections",
    subtitle = "Joint cluster assignments mapped back onto each section",
    caption  = "Module 16 · Spatial Transcriptomics")
ggsave(file.path(OUT_DIR, "Mod16_C8_integrated_spatial_clusters.png"), p_int_spatial, width = 10, height = 5, dpi = 300)

saveRDS(brain, file.path(OBJ_DIR, "brain_spatial_integrated.rds"))
cat("Wrote", file.path(OBJ_DIR, "brain_spatial_integrated.rds"), "with", ncol(brain), "spots\n")
cat("Wrote Mod16 figures/tables to", OUT_DIR, "\n")
