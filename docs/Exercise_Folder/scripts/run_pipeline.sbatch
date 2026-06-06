#!/bin/bash
#SBATCH --account=<myPIRG>
#SBATCH --partition=compute
#SBATCH --job-name=ifnb_pipeline
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G

# ---------------------------------------------------------------------------
# Example batch script: run the WHOLE core ifnb analysis pipeline (steps
# 01-08) end-to-end as a SINGLE SLURM job. Each script reads the .rds/.csv
# objects the previous one wrote, so they must run in order — this wrapper
# just runs them sequentially in one allocation, which is the simplest way
# to reproduce the entire laptop series (Modules 01-08) on Talapas.
#
# This is the script referenced by Module 09 (SLURM Basics) as the example
# of "how to run the series of .R scripts". For finer control — running one
# step at a time, or right-sizing memory per step — use run_rscript.sbatch
# instead (see Module 10, the Analysis Pipeline page).
#
# PREREQUISITES (one time, from this scripts/ folder):
#   mkdir -p logs
#   # pre-cache ifnb on a login node (compute nodes may lack internet):
#   Rscript -e 'invisible(muscData::Kang18_8vs8())'
#
# Submit:
#   sbatch run_pipeline.sbatch
#
# Watch / verify:
#   squeue -u $USER
#   seff <jobid>          # confirm State: COMPLETED (not OUT_OF_MEMORY/TIMEOUT)
#   tail -f logs/ifnb_pipeline_*.out
# ---------------------------------------------------------------------------

set -euo pipefail

module purge
module load R/4.4.0 || true   # adjust to a module available on Talapas

# The ordered core pipeline. Each step's output is the next step's input.
STEPS=(
  01_qc_preprocessing.R
  02_dimreduction_clustering.R
  03_markers_annotation.R
  04_reference_annotation.R
  05_integration.R
  06_pseudobulk_de.R
  07_functional_analysis.R
  08_differential_abundance.R
)

for step in "${STEPS[@]}"; do
  echo "==== $(date '+%F %T')  RUNNING ${step} ===="
  Rscript "${step}"
  echo "==== $(date '+%F %T')  FINISHED ${step} ===="
done

echo "==== Pipeline complete: all 8 steps finished on $(hostname) ===="
