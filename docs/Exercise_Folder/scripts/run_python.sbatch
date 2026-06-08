#!/bin/bash
#SBATCH --account=<myPIRG>
#SBATCH --partition=compute
#SBATCH --job-name=talapas_python
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Generic SLURM wrapper for the PYTHON scripts in this folder — currently the
# scVelo RNA-velocity step, 14_rna_velocity_scvelo.py (the Python sibling of the
# R script 14_trajectory_cellcomm.R).
#
# PREREQUISITE: the logs/ directory must exist before you submit, or SLURM
# cannot open --output=logs/... and the job is rejected. From this folder:
#   mkdir -p logs
#
# Usage:
#   sbatch --job-name=velocity run_python.sbatch 14_rna_velocity_scvelo.py
# Pass the script's options via environment variables, e.g.:
#   LOOM=../data/ifnb_velocyto.loom OUT_DIR=../output/Mod14 \
#     sbatch --job-name=velocity run_python.sbatch 14_rna_velocity_scvelo.py

set -euo pipefail

SCRIPT=${1:?usage: run_python.sbatch <script.py>}

module purge
module load miniconda3 || module load python || true  # adjust to a module on Talapas
# Activate the environment that has scvelo + scanpy installed (adjust the name):
# conda activate scvelo   ||   source activate scvelo

python "${SCRIPT}"
