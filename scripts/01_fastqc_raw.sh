#!/usr/bin/env bash
# ==============================================================================
# Script: 01_fastqc_raw.sh
# Description: Quality assessment of raw sequencing reads using FastQC & MultiQC
# ==============================================================================

set -euo pipefail

RAW_DIR="00_raw_fastq"
OUTDIR="01_fastqc_raw"
THREADS=8

mkdir -p "${OUTDIR}"

echo "=========================================================="
echo "Phase 1: Running FastQC on Raw FASTQ Files"
echo "=========================================================="

fastqc -t "${THREADS}" \
    -o "${OUTDIR}" \
    "${RAW_DIR}"/*.fastq*

# Compile MultiQC summary report if multiqc is installed
if command -v multiqc &>/dev/null; then
    echo "Running MultiQC on raw FastQC reports..."
    multiqc "${OUTDIR}" -o "${OUTDIR}/multiqc_report"
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Raw FastQC analysis finished."
