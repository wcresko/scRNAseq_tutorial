# Code checks

## `check_code.R` — static checker (no packages, no data, no network)

Run from the repo root:

```bash
Rscript tools/check_code.R
```

It validates, across all 14 `Exercise_Folder/Tutorial_*.qmd` and all
`Exercise_Folder/Talapas/scripts/*.R`:

1. **Syntax** — `parse()` of every `.R` file and every `{r}` chunk (unbalanced
   parens, stray commas, malformed calls).
2. **Library coverage** — flags a script that calls `ggsave`/`labs`/`write_csv`/
   `plot_annotation`/… without loading the package that provides it
   (the `could not find function "…"` class of runtime error).
3. **Pipeline hand-offs** — every `readRDS()`/`read_csv()` of an `OBJ_DIR`/`data`
   file must be produced by an *earlier* pipeline script
   (the `cannot open the connection` class of runtime error).

Exit code is `0` on PASS, `1` if any issue is found, so it is CI-friendly.

### What it does NOT catch

It does **not execute** the analysis, so it cannot surface package-API changes
(e.g. a Seurat v5 argument rename) or data-dependent runtime errors. Those need
an environment with the actual packages installed.

## Getting a true runtime check

The Claude Code on-the-web container used to author these files cannot reach
CRAN or Bioconductor (the network policy returns HTTP 403), so Seurat / DESeq2 /
Signac / … cannot be installed here and the pipelines cannot be executed.

To run the code for real, pick whichever is convenient:

- **Laptop / RStudio** — install the packages listed in each tutorial's setup
  chunk, then knit the `.qmd` with `eval: true` (or `source()` the `.R`).
- **Talapas (HPC)** — the `Talapas/scripts/*.R` are written for this; submit via
  `run_rscript.sbatch` (they download `ifnb` once through `ExperimentHub`).
- **Open the environment's network policy** — re-create the web environment with
  a policy that allows CRAN/Bioconductor, then `Rscript` can install and run.
  See https://code.claude.com/docs/en/claude-code-on-the-web (network policy).
- **GitHub Actions** — a hosted runner has open network and can install the
  packages + render the tutorials to catch runtime errors in CI.
