#!/usr/bin/env bash
#SBATCH --job-name=shotgun_metagenome_SRR24315725
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err

# ==============================================================================
# SLURM Job Submission Script for High-Throughput Shotgun Metagenomics
# Tested for Zaratan / HPC Clusters with Slurm Workload Manager
# ==============================================================================

set -euo pipefail

echo "Job ID: $SLURM_JOB_ID"
echo "Running on Host: $(hostname)"
echo "Started at: $(date)"

# Load Environment Modules (Uncomment/adjust for your cluster)
# module purge
# module load fastqc/0.12.1
# module load fastp/0.23.4
# module load bowtie2/2.5.2
# module load samtools/1.19
# module load kraken2/2.1.3
# module load bracken/2.9
# module load humann/3.8
# module load R/4.3.3

# Or activate Conda / Mamba Environment:
# source ~/.bashrc
# conda activate metagenomics_env

export THREADS=$SLURM_CPUS_PER_TASK

# Set database paths on HPC scratch / shared storage:
export BOWTIE2_HUMAN_INDEX="/shared/databases/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as"
export KRAKEN2_DB="/shared/databases/kraken2/k2_standard_2024"
export HUMANN_NUC_DB="/shared/databases/humann/chocophlan"
export HUMANN_PROT_DB="/shared/databases/humann/uniref50"

# Execute master pipeline
chmod +x ./run_metagenomics_pipeline.sh
./run_metagenomics_pipeline.sh

echo "Finished at: $(date)"
