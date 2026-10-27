#!/usr/bin/env bash
#SBATCH --job-name=shotgun_crc_vs_healthy
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL

# ==============================================================================
# SLURM Batch Job Script for Shotgun Metagenomics Pipeline
# Optimized for HPC Clusters (e.g. UMD Zaratan, NIH Biowulf, Academic Clusters)
# ==============================================================================

set -euo pipefail

mkdir -p logs

echo "Job started on $(hostname) at $(date)"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK, Memory: 64G"

# Load modules or activate conda environment
# Example for Zaratan / standard cluster:
# module load fastqc fastp bowtie2/2.5.1 kraken2/2.1.2 R/4.3.0
# Or activate conda environment:
source "$(conda info --base)/etc/profile.d/conda.sh" || true
conda activate shotgun_env || true

# Export database paths (adjust to your cluster scratch/reference directory)
export BOWTIE2_INDEX_DIR="/scratch/zt1/project/references/bowtie2/GRCh38"
export KRAKEN2_DB="/scratch/zt1/project/references/kraken2/k2_standard"
export CHOCOPHLAN_DIR="/scratch/zt1/project/references/humann3/chocophlan"
export UNIREF_DIR="/scratch/zt1/project/references/humann3/uniref"

# Execute master workflow
bash scripts/run_pipeline.sh

echo "Job finished at $(date)"
