#!/usr/bin/env bash
# ==============================================================================
# Script: 02_fastp_trim.sh
# Description: Quality filtering, adapter trimming, and polyG tail removal via fastp
# ==============================================================================

set -euo pipefail

IN_DIR="00_raw_fastq"
OUT_DIR="02_fastp_trimmed"
QC_POST_DIR="02_fastqc_post_trim"
THREADS=8

mkdir -p "${OUT_DIR}" "${QC_POST_DIR}"

SAMPLES=(
    "SRR24315726"
    "SRR24315725"
    "SRR24315724"
    "SRR24315723"
)

echo "=========================================================="
echo "Phase 2: Running fastp Trimming and Filtering"
echo "=========================================================="

for ACC in "${SAMPLES[@]}"; do
    # Handle both .fastq and .fastq.gz
    if [[ -f "${IN_DIR}/${ACC}_1.fastq.gz" ]]; then
        R1="${IN_DIR}/${ACC}_1.fastq.gz"
        R2="${IN_DIR}/${ACC}_2.fastq.gz"
    elif [[ -f "${IN_DIR}/${ACC}_1.fastq" ]]; then
        R1="${IN_DIR}/${ACC}_1.fastq"
        R2="${IN_DIR}/${ACC}_2.fastq"
    elif [[ -f "${ACC}_1.fastq" ]]; then
        R1="${ACC}_1.fastq"
        R2="${ACC}_2.fastq"
    else
        echo "Error: FASTQ files for ${ACC} not found!" >&2
        continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Trimming sample ${ACC}..."

    fastp \
        --in1 "${R1}" \
        --in2 "${R2}" \
        --out1 "${OUT_DIR}/${ACC}_1.clean.fastq.gz" \
        --out2 "${OUT_DIR}/${ACC}_2.clean.fastq.gz" \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --cut_mean_quality 20 \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 30 \
        --n_base_limit 5 \
        --length_required 50 \
        --trim_poly_g \
        --thread "${THREADS}" \
        --html "${OUT_DIR}/${ACC}_fastp.html" \
        --json "${OUT_DIR}/${ACC}_fastp.json"

    # Post-trimming FastQC check
    fastqc -t "${THREADS}" \
        -o "${QC_POST_DIR}" \
        "${OUT_DIR}/${ACC}_1.clean.fastq.gz" \
        "${OUT_DIR}/${ACC}_2.clean.fastq.gz"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] fastp trimming and post-QC completed."
