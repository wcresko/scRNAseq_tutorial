#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# build_workshop_rds.R
#
# Builds the two core-module datasets as plain .rds files so students who
# cannot get the Bioconductor data packages installed (muscData / airway)
# can skip the package step entirely: they download one .rds and read it in
# with readRDS() as their first step.
#
#   ifnb_raw.rds  - the assembled `ifnb` Seurat object, exactly as it leaves
#                   Step 1 of Tutorial 01 (CTRL + STIM, cleaned gene symbols,
#                   `stim` and `seurat_annotations` metadata) but BEFORE any
#                   QC. Drop-in replacement for the Kang18_8vs8() block.
#   airway_raw.rds    - the `airway` RangedSummarizedExperiment, identical to
#                   `data(airway)`. Drop-in replacement for that line in
#                   Tutorial 06.
#
# This is meant to be run by the instructor (it needs muscData + airway +
# internet once). render_all.sh runs it in the background; it can also be
# run on its own:
#
#     Rscript tools/build_workshop_rds.R [OUTPUT_DIR]
#
# OUTPUT_DIR defaults to docs/rds (so the files publish with the site and
# are downloadable from the Datasets page). Failures here are non-fatal to
# the website build: the page also documents the package-based path.
# -------------------------------------------------------------------------

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "docs/rds"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("=== build_workshop_rds.R ===")
message("Writing .rds datasets to: ", normalizePath(out_dir, mustWork = FALSE))

# GitHub blocks pushes of files larger than 100 MB. Warn the instructor so a
# too-large object can be re-hosted (Release asset / cloud drive) instead.
GITHUB_FILE_LIMIT_MB <- 95

report_size <- function(path) {
  if (!file.exists(path)) {
    message("  !! expected output missing: ", path)
    return(invisible())
  }
  size_mb <- file.info(path)$size / 1024^2
  message(sprintf("  wrote %s (%.1f MB)", basename(path), size_mb))
  if (size_mb > GITHUB_FILE_LIMIT_MB) {
    message(sprintf(
      "  !! %s is %.1f MB (> %d MB). GitHub will reject the push. Host it on a\n     GitHub Release or cloud drive instead and update the link on Datasets.qmd.",
      basename(path), size_mb, GITHUB_FILE_LIMIT_MB))
  }
}

# -------------------------------------------------------------------------
# 1. ifnb -> ifnb_raw.rds
#    Mirrors Tutorial_01_QC_Preprocessing.qmd Step 1 (label: M1-load) exactly
#    so the saved object is byte-for-byte what students would build in class.
# -------------------------------------------------------------------------
build_ifnb <- function() {
  suppressPackageStartupMessages({
    library(Seurat)
    library(muscData)
    library(SingleCellExperiment)
  })

  # Fetch the Kang et al. 2017 PBMC dataset (CTRL + IFN-beta STIM) from
  # Bioconductor's ExperimentHub. First call downloads ~25 MB and caches it.
  sce <- Kang18_8vs8()

  # Keep only true singlets that have an author-assigned cell-type label.
  sce <- sce[, sce$multiplets == "singlet" & !is.na(sce$cell)]

  # Gene-name cleanup: muscData stores each feature as SYMBOL_ENSEMBLID. Keep
  # only the readable SYMBOL (strip a trailing _ENSEMBLID), preserving real
  # dashed symbols (HLA-A, MT-CO1) and de-duplicating with make.unique().
  cm    <- counts(sce)
  orig  <- rownames(cm)
  sym   <- sub("_ENS[A-Z]*[0-9]+(\\.[0-9]+)?$", "", orig)
  blank <- is.na(sym) | sym == "" | grepl("^ENS[A-Z]*[0-9]+", sym)
  sym[blank] <- orig[blank]
  rownames(cm) <- make.unique(gsub("_", "-", sym))

  # Assemble the Seurat object with the metadata column names the rest of the
  # tutorial series expects.
  ifnb <- CreateSeuratObject(
    counts    = cm,
    meta.data = as.data.frame(colData(sce))
  )
  ifnb$stim <- factor(toupper(ifnb$stim), levels = c("CTRL", "STIM"))
  ifnb$seurat_annotations <- factor(ifnb$cell)

  path <- file.path(out_dir, "ifnb_raw.rds")
  saveRDS(ifnb, file = path, compress = "xz")   # xz for the smallest file
  report_size(path)
}

# -------------------------------------------------------------------------
# 2. airway -> airway_raw.rds
#    Identical to `library(airway); data(airway)` in Tutorial 06 Part A.
# -------------------------------------------------------------------------
build_airway <- function() {
  suppressPackageStartupMessages(library(airway))
  data(airway)
  path <- file.path(out_dir, "airway_raw.rds")
  saveRDS(airway, file = path, compress = "xz")
  report_size(path)
}

ok_ifnb <- tryCatch({ build_ifnb(); TRUE },
  error = function(e) { message("  !! ifnb build failed: ", conditionMessage(e)); FALSE })

ok_airway <- tryCatch({ build_airway(); TRUE },
  error = function(e) { message("  !! airway build failed: ", conditionMessage(e)); FALSE })

message("=== build_workshop_rds.R done (ifnb: ", if (ok_ifnb) "ok" else "FAILED",
        ", airway: ", if (ok_airway) "ok" else "FAILED", ") ===")

# Exit non-zero only if BOTH failed, so render_all can surface a real problem
# without breaking the site build for a transient single-dataset hiccup.
if (!ok_ifnb && !ok_airway) quit(status = 1)
