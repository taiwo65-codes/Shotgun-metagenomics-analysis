#!/usr/bin/env bash
# ==============================================================================
# Script: 05_humann3_profiling.sh
# Description: Functional metabolic and pathway profiling using HUMAnN3
# ==============================================================================

set -euo pipefail

IN_DIR="03_decontaminated_reads"
HUMANN_OUT="05_humann3_output"
JOINED_DIR="05_humann3_joined"
THREADS=16

# HUMAnN3 Database Paths
CHOCOPHLAN_DB="${CHOCOPHLAN_DIR:-/databases/humann3/chocophlan}"
UNIREF_DB="${UNIREF_DIR:-/databases/humann3/uniref}"

mkdir -p "${HUMANN_OUT}" "${JOINED_DIR}"

SAMPLES=(
    "SRR24315726"
    "SRR24315725"
    "SRR24315724"
    "SRR24315723"
)

echo "=========================================================="
echo "Phase 5: Functional Profiling with HUMAnN3"
echo "ChocoPhlAn: ${CHOCOPHLAN_DB}"
echo "UniRef DB:  ${UNIREF_DB}"
echo "=========================================================="

for ACC in "${SAMPLES[@]}"; do
    R1="${IN_DIR}/${ACC}_1.decontam.fastq.gz"
    R2="${IN_DIR}/${ACC}_2.decontam.fastq.gz"

    if [[ ! -f "${R1}" ]]; then
        R1="02_fastp_trimmed/${ACC}_1.clean.fastq.gz"
        R2="02_fastp_trimmed/${ACC}_2.clean.fastq.gz"
    fi

    # Concatenate paired-end reads for HUMAnN input
    MERGED_FQ="${HUMANN_OUT}/${ACC}_concat.fastq.gz"
    if [[ ! -f "${MERGED_FQ}" ]]; then
        echo "Concatenating pairs for ${ACC}..."
        cat "${R1}" "${R2}" > "${MERGED_FQ}"
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running HUMAnN3 on ${ACC}..."
    humann \
        --input "${MERGED_FQ}" \
        --output "${HUMANN_OUT}/${ACC}" \
        --threads "${THREADS}" \
        --nucleotide-database "${CHOCOPHLAN_DB}" \
        --protein-database "${UNIREF_DB}" \
        --remove-temp-output

    # Normalize gene families and pathways to Copies Per Million (CPM)
    echo "Normalizing abundances for ${ACC} to CPM..."
    humann_renorm_table \
        --input "${HUMANN_OUT}/${ACC}/${ACC}_concat_pathabundance.tsv" \
        --output "${HUMANN_OUT}/${ACC}/${ACC}_pathabundance_cpm.tsv" \
        --units cpm

    humann_renorm_table \
        --input "${HUMANN_OUT}/${ACC}/${ACC}_concat_genefamilies.tsv" \
        --output "${HUMANN_OUT}/${ACC}/${ACC}_genefamilies_cpm.tsv" \
        --units cpm
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Joining sample tables across cohort..."
# Join Pathway Abundance tables
humann_join_tables \
    --input "${HUMANN_OUT}" \
    --output "${JOINED_DIR}/cohort_pathabundance_cpm.tsv" \
    --file_name "pathabundance_cpm"

# Join Gene Family tables
humann_join_tables \
    --input "${HUMANN_OUT}" \
    --output "${JOINED_DIR}/cohort_genefamilies_cpm.tsv" \
    --file_name "genefamilies_cpm"

# Split stratified vs unstratified pathway tables
humann_split_stratified_table \
    --input "${JOINED_DIR}/cohort_pathabundance_cpm.tsv" \
    --output "${JOINED_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] HUMAnN3 functional profiling complete."
