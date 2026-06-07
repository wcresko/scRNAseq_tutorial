#!/bin/bash
#SBATCH --account=<myPIRG>
#SBATCH --partition=compute
#SBATCH --job-name=cellranger_atac
#SBATCH --output=logs/cellranger_atac_%j.out
#SBATCH --error=logs/cellranger_atac_%j.err
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G

# Module 12 — Raw scATAC-seq processing with cellranger-atac count (Talapas)
# Edit <PIRG>, <user>, <SAMPLE>, and paths below to match your environment.

set -euo pipefail

PIRG=<myPIRG>
USER=<your-duckid>
PROJ=/projects/${PIRG}/${USER}/atac
REF=${PROJ}/refs/refdata-cellranger-arc-GRCh38-2024-A
FASTQS=${PROJ}/raw_fastqs
RUNS=${PROJ}/cellranger_atac_runs
SAMPLE=atac_pbmc_1k_v1   # change to your sample ID

mkdir -p ${PROJ} ${FASTQS} ${RUNS} ${PROJ}/refs logs

# ATAC reference (note: 'arc' build, NOT the RNA 'gex' reference): download once per group.
if [ ! -d "${REF}" ]; then
  cd ${PROJ}/refs
  wget -q https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-GRCh38-2024-A.tar.gz
  tar -xzf refdata-cellranger-arc-GRCh38-2024-A.tar.gz
fi

# Raw scATAC FASTQs: download a small public test set (PBMC 1k v1 ATAC, ~3 GB).
if [ ! -d "${FASTQS}/atac_pbmc_1k_v1_fastqs" ]; then
  cd ${FASTQS}
  wget -q https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-atac/1.2.0/atac_pbmc_1k_v1/atac_pbmc_1k_v1_fastqs.tar
  tar -xf atac_pbmc_1k_v1_fastqs.tar
fi

# cellranger-atac: load module or set PATH.
module load cellranger-atac/2.2.0 || \
  export PATH=${PROJ}/software/cellranger-atac-2.2.0:${PATH}

# Run cellranger-atac count.
cd ${RUNS}
cellranger-atac count \
  --id=${SAMPLE} \
  --reference=${REF} \
  --fastqs=${FASTQS}/atac_pbmc_1k_v1_fastqs \
  --sample=${SAMPLE} \
  --localcores=16 \
  --localmem=120

echo "Done. fragments.tsv.gz + peaks at ${RUNS}/${SAMPLE}/outs/ (feed Module 15 / Signac)"
