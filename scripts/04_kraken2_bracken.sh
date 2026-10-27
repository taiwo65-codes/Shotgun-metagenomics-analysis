#!/usr/bin/env bash
# ==============================================================================
# Script: 04_kraken2_bracken.sh
# Description: Taxonomic profiling using Kraken2 followed by Bracken abundance estimation
# ==============================================================================

set -euo pipefail

IN_DIR="03_decontaminated_reads"
KRAKEN_DIR="04_kraken2_output"
BRACKEN_DIR="04_bracken_output"
THREADS=16

# Kraken2 Database Path (Standard DB, PlusPF, or MiniKraken)
# e.g., /databases/kraken2/k2_standard_20240904 or /scratch/zaratandb/kraken2
KRAKEN_DB="${KRAKEN2_DB:-/databases/kraken2/standard}"
READ_LEN=150
TAX_LEVEL="S"   # S for Species, G for Genus, P for Phylum
THRESHOLD=10    # Minimum reads required to re-estimate abundance

mkdir -p "${KRAKEN_DIR}" "${BRACKEN_DIR}"

SAMPLES=(
    "SRR24315726"
    "SRR24315725"
    "SRR24315724"
    "SRR24315723"
)

echo "=========================================================="
echo "Phase 4: Taxonomic Profiling with Kraken2 & Bracken"
echo "Database: ${KRAKEN_DB}"
echo "Level: ${TAX_LEVEL} (Species), Read length: ${READ_LEN}"
echo "=========================================================="

for ACC in "${SAMPLES[@]}"; do
    R1="${IN_DIR}/${ACC}_1.decontam.fastq.gz"
    R2="${IN_DIR}/${ACC}_2.decontam.fastq.gz"

    # Fallback to trimmed if decontaminated not present
    if [[ ! -f "${R1}" ]]; then
        R1="02_fastp_trimmed/${ACC}_1.clean.fastq.gz"
        R2="02_fastp_trimmed/${ACC}_2.clean.fastq.gz"
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Kraken2 on ${ACC}..."
    kraken2 \
        --db "${KRAKEN_DB}" \
        --threads "${THREADS}" \
        --paired \
        --gzip-compressed \
        --report "${KRAKEN_DIR}/${ACC}.kraken2.report" \
        --output "${KRAKEN_DIR}/${ACC}.kraken2.output" \
        "${R1}" "${R2}"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Bracken Bayesian abundance re-estimation on ${ACC}..."
    bracken \
        -d "${KRAKEN_DB}" \
        -i "${KRAKEN_DIR}/${ACC}.kraken2.report" \
        -o "${BRACKEN_DIR}/${ACC}.bracken.species.txt" \
        -w "${BRACKEN_DIR}/${ACC}.bracken.report" \
        -r "${READ_LEN}" \
        -l "${TAX_LEVEL}" \
        -t "${THRESHOLD}"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Merging Bracken species outputs into abundance matrix..."
combine_bracken_outputs \
    --files "${BRACKEN_DIR}"/*.bracken.species.txt \
    --output "${BRACKEN_DIR}/merged_bracken_species_abundance.txt"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Taxonomic profiling completed successfully."
