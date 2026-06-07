#!/bin/bash
#SBATCH --account=<myPIRG>
#SBATCH --partition=compute
#SBATCH --job-name=cellranger_count
#SBATCH --output=logs/cellranger_%j.out
#SBATCH --error=logs/cellranger_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G

# Module 11 — Raw data download + cellranger count (Talapas)
# Edit <PIRG>, <user>, <SAMPLE>, and paths below to match your environment.

set -euo pipefail

PIRG=<myPIRG>
USER=<your-duckid>
PROJ=/projects/${PIRG}/${USER}/nsclc
REF=${PROJ}/refs/refdata-gex-GRCh38-2024-A
FASTQS=${PROJ}/raw_fastqs
RUNS=${PROJ}/cellranger_runs
SAMPLE=pbmc_1k_v3   # change to your sample ID

mkdir -p ${PROJ} ${FASTQS} ${RUNS} ${PROJ}/refs logs

# Reference: download once per group, share thereafter.
if [ ! -d "${REF}" ]; then
  cd ${PROJ}/refs
  wget -q https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
  tar -xzf refdata-gex-GRCh38-2024-A.tar.gz
fi

# Raw FASTQs: download a small public test set (PBMC 1k v3, ~5 GB).
if [ ! -d "${FASTQS}/pbmc_1k_v3_fastqs" ]; then
  cd ${FASTQS}
  wget -q https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
  tar -xf pbmc_1k_v3_fastqs.tar
fi

# Cell Ranger: load module or set PATH.
module load cellranger/9.0.0 || \
  export PATH=${PROJ}/software/cellranger-9.0.0:${PATH}

# Run cellranger count.
cd ${RUNS}
cellranger count \
  --id=${SAMPLE} \
  --transcriptome=${REF} \
  --fastqs=${FASTQS}/pbmc_1k_v3_fastqs \
  --sample=${SAMPLE} \
  --create-bam=true \
  --localcores=16 \
  --localmem=120

echo "Done. filtered_feature_bc_matrix at ${RUNS}/${SAMPLE}/outs/filtered_feature_bc_matrix"
