#!/bin/bash
#SBATCH --account=<myPIRG>
#SBATCH --partition=compute
#SBATCH --job-name=ifnb_bonus_pipeline
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G

# ---------------------------------------------------------------------------
# Example batch script: run the downstream BONUS analyses (Modules 13-17) as a
# single SLURM job — the parallel of run_core_pipeline.sbatch (Modules 01-08).
#
# IMPORTANT — unlike the core pipeline, these scripts are LARGELY STANDALONE:
# they do NOT chain on each other, and each needs its own input data:
#   13_wgcna.R            -> ../data/GSE152418_p20047_Study1_RawCounts.txt  (Module 13)
#   14_trajectory_cellcomm.R -> ../data/ifnb_annotated.rds   (from the core run; Module 14)
#   15_scatac_signac.R    -> ../data/atac_v1_pbmc_10k_*          (Module 15)
#   16_spatial.R          -> stxBrain via SeuratData            (Module 16)
#   17_fair_metadata.R    -> ../data/ifnb_annotated_final.rds (from the core run; Module 17)
# Because they are independent, this wrapper does NOT abort if one step fails
# (e.g. its dataset isn't downloaded) — it logs the failure and continues.
#
# RNA velocity (Module 14) is a PYTHON step (scVelo) and is NOT run here; submit
# it separately:  sbatch run_python.sbatch 14_rna_velocity_scvelo.py
# Raw-data Cell Ranger (Modules 11-12) are separate jobs: cellranger_RNA_count.sbatch
# and cellranger_atac_count.sbatch.
#
# PREREQUISITES (one time, from this scripts/ folder):
#   mkdir -p logs
#   # run the core pipeline first if you want 14 / 17 (they read its objects):
#   sbatch run_core_pipeline.sbatch
#   # and download the per-module datasets — see the Datasets page.
#
# Submit:
#   sbatch run_bonus_pipeline.sbatch
#
# Watch / verify:
#   squeue -u $USER
#   seff <jobid>
#   tail -f logs/ifnb_bonus_pipeline_*.out
# ---------------------------------------------------------------------------

set -uo pipefail   # NOTE: no -e — keep going if one standalone step fails

module purge
module load R/4.4.0 || true   # adjust to a module available on Talapas

# The bonus downstream R analyses (independent of one another).
STEPS=(
  13_wgcna.R
  14_trajectory_cellcomm.R
  15_scatac_signac.R
  16_spatial.R
  17_fair_metadata.R
)

fail=0
for step in "${STEPS[@]}"; do
  echo "==== $(date '+%F %T')  RUNNING ${step} ===="
  if Rscript "${step}"; then
    echo "==== $(date '+%F %T')  FINISHED ${step} ===="
  else
    echo "==== $(date '+%F %T')  ${step} FAILED (often a missing dataset) — continuing ===="
    fail=$((fail + 1))
  fi
done

echo "==== Bonus pipeline done on $(hostname): ${fail} of ${#STEPS[@]} step(s) failed ===="
echo "==== Reminder: RNA velocity is Python — sbatch run_python.sbatch 14_rna_velocity_scvelo.py ===="
