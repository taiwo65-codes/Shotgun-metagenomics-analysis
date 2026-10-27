#!/usr/bin/env bash
# ==============================================================================
# Script: 03_bowtie2_decontam.sh
# Description: Host DNA decontamination (Human GRCh38 removal) using Bowtie2
# ==============================================================================

set -euo pipefail

IN_DIR="02_fastp_trimmed"
OUT_DIR="03_decontaminated_reads"
LOG_DIR="03_bowtie2_logs"
THREADS=16

# Path to human GRCh38 Bowtie2 index
# If not yet downloaded, download pre-built index:
# wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip && unzip GRCh38_noalt_as.zip
INDEX_DIR="${BOWTIE2_INDEX_DIR:-/databases/bowtie2/GRCh38}"
INDEX_NAME="GRCh38_noalt_as"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

SAMPLES=(
    "SRR24315726"
    "SRR24315725"
    "SRR24315724"
    "SRR24315723"
)

echo "=========================================================="
echo "Phase 3: Human Host Decontamination via Bowtie2"
echo "Index: ${INDEX_DIR}/${INDEX_NAME}"
echo "=========================================================="

for ACC in "${SAMPLES[@]}"; do
    R1="${IN_DIR}/${ACC}_1.clean.fastq.gz"
    R2="${IN_DIR}/${ACC}_2.clean.fastq.gz"
    
    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "Clean reads for ${ACC} not found in ${IN_DIR}. Skipping."
        continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Aligning ${ACC} against human genome..."

    # Align against human genome; discard mapped reads and output unmapped pairs
    bowtie2 \
        -x "${INDEX_DIR}/${INDEX_NAME}" \
        -1 "${R1}" \
        -2 "${R2}" \
        --threads "${THREADS}" \
        --very-sensitive \
        --un-conc-gz "${OUT_DIR}/${ACC}_microbial_R%.fastq.gz" \
        -S /dev/null \
        2> "${LOG_DIR}/${ACC}_bowtie2_host_filter.log"

    # Rename to clean convention if needed
    if [[ -f "${OUT_DIR}/${ACC}_microbial_R1.fastq.gz" ]]; then
        mv "${OUT_DIR}/${ACC}_microbial_R1.fastq.gz" "${OUT_DIR}/${ACC}_1.decontam.fastq.gz"
        mv "${OUT_DIR}/${ACC}_microbial_R2.fastq.gz" "${OUT_DIR}/${ACC}_2.decontam.fastq.gz"
    fi

    # Parse mapping statistics
    TOTAL_PAIRS=$(grep "reads; of these:" "${LOG_DIR}/${ACC}_bowtie2_host_filter.log" | awk '{print $1}')
    ALIGN_RATE=$(grep "overall alignment rate" "${LOG_DIR}/${ACC}_bowtie2_host_filter.log" | awk '{print $1}')
    echo "Sample ${ACC}: Total pairs = ${TOTAL_PAIRS}, Host alignment rate = ${ALIGN_RATE}"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Host decontamination completed."
