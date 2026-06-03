#!/bin/bash
#SBATCH --account=<myPIRG>
#SBATCH --partition=compute
#SBATCH --job-name=slurm_demo
#SBATCH --output=logs/slurm_demo_%j.out
#SBATCH --error=logs/slurm_demo_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100M

# Minimal SLURM demo. Submit with: sbatch 01_slurm_demo.sbatch
# See the .qmd for what each #SBATCH directive does.
#
# PREREQUISITE: run `mkdir -p logs` once before submitting — SLURM will not
# create the logs/ directory that --output/--error point to.

echo "Hello from $(hostname) at $(date)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "SLURM_JOB_NAME=$SLURM_JOB_NAME"
sleep 10
echo "Done."
