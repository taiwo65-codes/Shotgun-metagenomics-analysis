#!/usr/bin/env bash
# ==============================================================================
# Script: 00_download_sra.sh
# Description: Download paired-end FASTQ data from NCBI SRA using fasterq-dump
# Author: Shotgun Metagenomics Pipeline
# ==============================================================================

set -euo pipefail

# Directory setup
OUTDIR="00_raw_fastq"
mkdir -p "${OUTDIR}"

# SRA Run Accessions:
# SRR24315726 -> SRX20110651 (CRC Patient 82)
# SRR24315725 -> SRX20110652 (CRC Patient 83)
# SRR24315724 -> SRX20110653 (Healthy Patient 55)
# SRR24315723 -> SRX20110654 (Healthy Patient 56)
SAMPLES=(
    "SRR24315726"
    "SRR24315725"
    "SRR24315724"
    "SRR24315723"
)

THREADS=8

echo "=========================================================="
echo "Starting SRA Download for Cohort (CRC vs Healthy)"
echo "Samples: ${SAMPLES[*]}"
echo "=========================================================="

for ACC in "${SAMPLES[@]}"; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloading ${ACC}..."
    
    if [[ -f "${OUTDIR}/${ACC}_1.fastq.gz" && -f "${OUTDIR}/${ACC}_2.fastq.gz" ]]; then
        echo "Files for ${ACC} already exist in ${OUTDIR}. Skipping."
        continue
    fi

    # Prefetch SRA file (optimal for large runs)
    prefetch "${ACC}" --output-directory "${OUTDIR}/tmp_sra"
    
    # Extract paired fastq files with fasterq-dump
    fasterq-dump "${OUTDIR}/tmp_sra/${ACC}/${ACC}.sra" \
        --split-files \
        --outdir "${OUTDIR}" \
        --threads "${THREADS}" \
        --temp "${OUTDIR}/tmp_sra"

    # Compress fastq files to save disk space
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Compressing ${ACC} FASTQs..."
    pigz -p "${THREADS}" "${OUTDIR}/${ACC}_1.fastq"
    pigz -p "${THREADS}" "${OUTDIR}/${ACC}_2.fastq"

    # Clean up temporary SRA files
    rm -rf "${OUTDIR}/tmp_sra/${ACC}"
done

rm -rf "${OUTDIR}/tmp_sra"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] SRA Download and compression completed successfully."
