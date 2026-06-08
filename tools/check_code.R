#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Static code checker for the workshop tutorials and Talapas pipeline scripts.
#
# Runs WITHOUT installing Seurat/Bioconductor or downloading any data, so it
# works anywhere base R is available. It catches the bug classes that fail
# before (or independently of) package + data loading:
#
#   1. SYNTAX      — parse() every .R script and every {r} chunk in the .qmd
#                    tutorials. Catches unbalanced parens, bad commas, typos.
#   2. LIBRARY     — flag scripts that call ggplot2/dplyr/readr/patchwork/...
#                    functions without loading the package that provides them
#                    ("could not find function" at runtime).
#   3. HAND-OFFS   — every readRDS()/read_csv() of a DATA_DIR/ data file must be
#                    produced by an earlier pipeline script's saveRDS/write_csv
#                    ("cannot open the connection" at runtime).
#
# It does NOT execute the analysis, so it cannot catch package-API or data-
# dependent runtime errors — for those, run in an environment that has the
# packages (laptop / Talapas), or see tools/README_checks.md.
#
# Usage (from repo root):  Rscript tools/check_code.R
# ---------------------------------------------------------------------------

root <- normalizePath(file.path(dirname(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
if (is.na(root) || !dir.exists(file.path(root, "Exercise_Folder"))) root <- getwd()
setwd(root)

SCRIPTS_DIR <- "Exercise_Folder/scripts"
fails <- 0L

## ----------------------------------------------------------------- helpers
strip_comments <- function(txt)
  paste(sub("#.*$", "", strsplit(txt, "\n", fixed = TRUE)[[1]]), collapse = "\n")

extract_chunks <- function(path) {
  lines <- readLines(path, warn = FALSE); chunks <- list()
  inchunk <- FALSE; buf <- character(); start <- 0L; lbl <- ""
  for (i in seq_along(lines)) {
    ln <- lines[i]
    if (!inchunk && grepl("^```\\{r", ln)) { inchunk <- TRUE; buf <- character(); start <- i; lbl <- ""; next }
    if (inchunk && grepl("^```\\s*$", ln)) {
      chunks[[length(chunks) + 1]] <- list(code = paste(buf, collapse = "\n"), start = start, label = lbl)
      inchunk <- FALSE; next
    }
    if (inchunk) {
      m <- regmatches(ln, regexpr("(?<=#\\| label:\\s).*", ln, perl = TRUE))
      if (length(m)) lbl <- trimws(m)
      buf <- c(buf, ln)
    }
  }
  chunks
}

try_parse <- function(code, label) {
  err <- tryCatch({ parse(text = code); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(err)) {
    cat(sprintf("  [PARSE ERROR] %s\n      %s\n", label, gsub("\n", "\n      ", trimws(err))))
    return(FALSE)
  }
  TRUE
}

## ----------------------------------------------------------------- 1. SYNTAX
cat("==== 1. SYNTAX ====\n")
rfiles <- list.files(SCRIPTS_DIR, pattern = "\\.R$", full.names = TRUE)
for (f in rfiles)
  if (!try_parse(paste(readLines(f, warn = FALSE), collapse = "\n"), f)) fails <- fails + 1L
nchunks <- 0L
for (f in list.files("Exercise_Folder", "Tutorial_.*\\.qmd$", full.names = TRUE))
  for (c in extract_chunks(f)) {
    nchunks <- nchunks + 1L
    lab <- sprintf("%s : chunk @L%d (label: %s)", basename(f), c$start,
                   if (nzchar(c$label)) c$label else "<none>")
    if (!try_parse(c$code, lab)) fails <- fails + 1L
  }
cat(sprintf("   parsed %d .R files + %d .qmd chunks\n", length(rfiles), nchunks))

## ----------------------------------------------------------------- 2. LIBRARY
cat("\n==== 2. LIBRARY COVERAGE ====\n")
NEED <- list(
  ggplot2   = "\\b(ggsave|ggtitle|labs|theme|element_text|element_blank|aes|geom_[a-z_]+|scale_[a-z_]+|guides|guide_legend|xlim|ylim|ggplot)\\b",
  patchwork = "\\b(plot_annotation|wrap_plots|wrap_elements|plot_layout|plot_spacer)\\b",
  readr     = "\\b(write_csv|read_csv|write_tsv|read_tsv)\\b",
  dplyr     = "\\b(mutate|filter|group_by|summarise|arrange|distinct|slice_max|slice_min|inner_join|bind_rows|pull|ungroup|across|transmute)\\b",
  tibble    = "\\b(tibble|enframe|as_tibble|tribble|column_to_rownames|rownames_to_column)\\b",
  purrr     = "\\b(map_dfr|map_chr|map_dbl|set_names|compact|imap)\\b",
  tidyr     = "\\b(pivot_longer|pivot_wider|unnest)\\b")
ATTACH <- list(tidyverse = c("ggplot2","dplyr","tibble","readr","purrr","tidyr","stringr","forcats"))
loaded <- function(txt) {
  pk <- unique(regmatches(txt, gregexpr("(?<=library\\()[A-Za-z0-9._]+|(?<=require\\()[A-Za-z0-9._]+", txt, perl = TRUE))[[1]])
  for (b in names(ATTACH)) if (b %in% pk) pk <- union(pk, ATTACH[[b]])
  pk
}
libok <- TRUE
for (f in list.files(SCRIPTS_DIR, "^[0-9].*\\.R$", full.names = TRUE)) {
  raw <- paste(readLines(f, warn = FALSE), collapse = "\n"); code <- strip_comments(raw); pk <- loaded(raw)
  miss <- names(NEED)[vapply(names(NEED), function(p)
    grepl(NEED[[p]], code, perl = TRUE) && !(p %in% pk), logical(1))]
  if (length(miss)) { libok <- FALSE; fails <- fails + 1L
    cat(sprintf("  [MISSING LIB] %s -> %s\n", basename(f), paste(miss, collapse = ", "))) }
}
if (libok) cat("   all pipeline scripts load the packages they use\n")

## ----------------------------------------------------------------- 3. HAND-OFFS
cat("\n==== 3. PIPELINE HAND-OFFS ====\n")
EXTERNAL <- c("atac_v1_pbmc_10k_singlecell.csv")  # known external (10x) inputs
num <- function(p) as.integer(sub(".*/(\\d+)_.*", "\\1", p))
sc <- sort(list.files(SCRIPTS_DIR, "^[0-9].*\\.R$", full.names = TRUE))
# pull every quoted *.rds/*.csv/*.tsv token from lines that read or write files
toks <- function(line) {
  m <- regmatches(line, gregexpr('["\']([^"\']+\\.(?:rds|csv|tsv))["\']', line, perl = TRUE))[[1]]
  basename(gsub('["\']', "", m))
}
W <- list(); R <- list()
for (s in sc) {
  n <- as.character(num(s)); lines <- readLines(s, warn = FALSE)
  w <- character(); rd <- character()
  for (ln in lines) {
    if (grepl("saveRDS|write_csv|write_tsv|write\\.csv", ln)) w <- c(w, toks(ln))
    if (grepl("readRDS|read_csv|read_tsv|read\\.csv", ln))    rd <- c(rd, toks(ln))
  }
  W[[n]] <- unique(w); R[[n]] <- unique(rd)
}
firstwrite <- list()
for (n in sort(as.integer(names(W)))) for (f in W[[as.character(n)]]) if (is.null(firstwrite[[f]])) firstwrite[[f]] <- n
hoff <- TRUE
for (n in sort(as.integer(names(R)))) for (f in R[[as.character(n)]]) {
  if (f %in% EXTERNAL) next
  src <- firstwrite[[f]]
  if (is.null(src)) { cat(sprintf("  [UNRESOLVED READ] script %02d reads '%s' (no producer)\n", n, f)); hoff <- FALSE; fails <- fails + 1L }
  else if (src >= n) { cat(sprintf("  [ORDER BUG] script %02d reads '%s' written only by %02d\n", n, f, src)); hoff <- FALSE; fails <- fails + 1L }
}
if (hoff) cat("   every pipeline input resolves to an earlier producer\n")

## ----------------------------------------------------------------- summary
cat(sprintf("\n==== RESULT: %s (%d issue%s) ====\n",
            if (fails == 0) "PASS" else "ISSUES FOUND", fails, if (fails == 1) "" else "s"))
quit(status = if (fails == 0) 0 else 1)
