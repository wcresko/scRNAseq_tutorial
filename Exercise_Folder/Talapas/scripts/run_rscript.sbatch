#!/bin/bash
#SBATCH --account=<myPIRG>
#SBATCH --partition=compute
#SBATCH --job-name=talapas_rscript
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Generic SLURM wrapper for the R scripts in this folder.
#
# PREREQUISITE: the logs/ directory must exist before you submit, or SLURM
# cannot open --output=logs/... and the job is rejected. From this folder:
#   mkdir -p logs
#
# Usage:
#   sbatch --job-name=qc run_rscript.sbatch 01_qc_preprocessing.R
#   sbatch --job-name=de run_rscript.sbatch 06_pseudobulk_de.R

set -euo pipefail

SCRIPT=${1:?usage: run_rscript.sbatch <script.R>}

module purge
module load R/4.4.0 || true   # adjust to a module available on Talapas

Rscript "${SCRIPT}"
