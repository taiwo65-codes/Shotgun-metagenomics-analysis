#!/usr/bin/env bash
# ==============================================================================
# Shotgun Metagenomics Analysis Pipeline for Colorectal Cancer (CRC) Cohorts
# Sample Accession: SRR24315725 (and scalable to multi-sample cohorts)
# Workflow Steps:
#   1. FastQC (Raw quality control)
#   2. fastp (Adapter, polyG, and quality trimming)
#   3. FastQC (Clean quality control)
#   4. Bowtie2 (Human host DNA decontamination against GRCh38)
#   5. Kraken2 & Bracken (Taxonomic profiling & Bayesian abundance re-estimation)
#   6. HUMAnN 3 (Functional profiling: Gene families, MetaCyc pathways)
#   7. MaAsLin 3 (Multivariable association analysis with clinical covariates)
# ==============================================================================

set -euo pipefail

# ----------------- CONFIGURATION & PATHS -----------------
SAMPLE="SRR24315725"
THREADS=8

# Working directories
RAW_DIR="."
OUT_DIR="./results"
FASTQC_RAW_DIR="${OUT_DIR}/01_fastqc_raw"
FASTP_DIR="${OUT_DIR}/02_fastp_trimmed"
FASTQC_TRIM_DIR="${OUT_DIR}/03_fastqc_trimmed"
DECONTAM_DIR="${OUT_DIR}/04_host_removal_bowtie2"
TAXA_DIR="${OUT_DIR}/05_taxonomic_kraken2_bracken"
FUNC_DIR="${OUT_DIR}/06_functional_humann3"
MAASLIN_DIR="${OUT_DIR}/07_maaslin3"

# Reference Database Paths (Set to your local or HPC paths)
BOWTIE2_HUMAN_INDEX="${BOWTIE2_HUMAN_INDEX:-/shared/databases/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as}"
KRAKEN2_DB="${KRAKEN2_DB:-/shared/databases/kraken2/standard}"
HUMANN_NUC_DB="${HUMANN_NUC_DB:-/shared/databases/humann/chocophlan}"
HUMANN_PROT_DB="${HUMANN_PROT_DB:-/shared/databases/humann/uniref50}"

mkdir -p "${FASTQC_RAW_DIR}" "${FASTP_DIR}" "${FASTQC_TRIM_DIR}" \
         "${DECONTAM_DIR}" "${TAXA_DIR}" "${FUNC_DIR}" "${MAASLIN_DIR}"

echo "======================================================================"
echo "Starting Shotgun Metagenomics Pipeline for Sample: ${SAMPLE}"
echo "Threads: ${THREADS}"
echo "======================================================================"

# ==============================================================================
# STEP 1: FastQC on Raw Reads
# ==============================================================================
echo "[STEP 1] Running FastQC on raw reads..."
fastqc -t "${THREADS}" \
  -o "${FASTQC_RAW_DIR}" \
  "${RAW_DIR}/${SAMPLE}_1.fastq" \
  "${RAW_DIR}/${SAMPLE}_2.fastq"

echo "[STEP 1] FastQC raw complete. Reports saved to ${FASTQC_RAW_DIR}"

# ==============================================================================
# STEP 2: fastp Trimming & Filtering
# ==============================================================================
echo "[STEP 2] Running fastp adapter trimming and quality filtering..."
fastp \
  -i "${RAW_DIR}/${SAMPLE}_1.fastq" \
  -I "${RAW_DIR}/${SAMPLE}_2.fastq" \
  -o "${FASTP_DIR}/${SAMPLE}_1.clean.fastq.gz" \
  -O "${FASTP_DIR}/${SAMPLE}_2.clean.fastq.gz" \
  --detect_adapter_for_pe \
  --cut_front \
  --cut_tail \
  --cut_mean_quality 20 \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 30 \
  --length_required 50 \
  --thread "${THREADS}" \
  --html "${FASTP_DIR}/${SAMPLE}_fastp.html" \
  --json "${FASTP_DIR}/${SAMPLE}_fastp.json"

echo "[STEP 2] fastp complete. Clean reads saved to ${FASTP_DIR}"

# ==============================================================================
# STEP 3: FastQC on Clean Reads
# ==============================================================================
echo "[STEP 3] Running FastQC on clean trimmed reads..."
fastqc -t "${THREADS}" \
  -o "${FASTQC_TRIM_DIR}" \
  "${FASTP_DIR}/${SAMPLE}_1.clean.fastq.gz" \
  "${FASTP_DIR}/${SAMPLE}_2.clean.fastq.gz"

echo "[STEP 3] Post-trimming FastQC complete."

# ==============================================================================
# STEP 4: Bowtie2 Alignment & Host Decontamination (GRCh38 Human Host Removal)
# ==============================================================================
echo "[STEP 4] Aligning against Human Genome (GRCh38) to remove host DNA..."
# In clinical stool/biopsy metagenomics, human DNA must be purged before microbial profiling.
# --un-conc-gz outputs the unmapped (microbial) paired reads directly.
bowtie2 \
  -p "${THREADS}" \
  -x "${BOWTIE2_HUMAN_INDEX}" \
  -1 "${FASTP_DIR}/${SAMPLE}_1.clean.fastq.gz" \
  -2 "${FASTP_DIR}/${SAMPLE}_2.clean.fastq.gz" \
  --very-sensitive \
  --un-conc-gz "${DECONTAM_DIR}/${SAMPLE}_nonhuman_%.fastq.gz" \
  --al-conc-gz "${DECONTAM_DIR}/${SAMPLE}_human_host_%.fastq.gz" \
  -S "${DECONTAM_DIR}/${SAMPLE}_human_align.sam"

# Compress / remove SAM to save disk space
samtools view -@ "${THREADS}" -bS "${DECONTAM_DIR}/${SAMPLE}_human_align.sam" > "${DECONTAM_DIR}/${SAMPLE}_human_align.bam"
rm -f "${DECONTAM_DIR}/${SAMPLE}_human_align.sam"

echo "[STEP 4] Host decontamination complete. Microbial reads saved to:"
echo "         ${DECONTAM_DIR}/${SAMPLE}_nonhuman_1.fastq.gz"
echo "         ${DECONTAM_DIR}/${SAMPLE}_nonhuman_2.fastq.gz"

# ==============================================================================
# STEP 5: Kraken2 & Bracken Taxonomic Profiling
# ==============================================================================
echo "[STEP 5a] Classifying reads with Kraken2..."
kraken2 \
  --threads "${THREADS}" \
  --db "${KRAKEN2_DB}" \
  --paired \
  --gzip-compressed \
  --report "${TAXA_DIR}/${SAMPLE}_kraken2_report.txt" \
  --output "${TAXA_DIR}/${SAMPLE}_kraken2_output.txt" \
  "${DECONTAM_DIR}/${SAMPLE}_nonhuman_1.fastq.gz" \
  "${DECONTAM_DIR}/${SAMPLE}_nonhuman_2.fastq.gz"

echo "[STEP 5b] Re-estimating species-level abundance with Bracken..."
bracken \
  -d "${KRAKEN2_DB}" \
  -i "${TAXA_DIR}/${SAMPLE}_kraken2_report.txt" \
  -o "${TAXA_DIR}/${SAMPLE}_bracken_species.txt" \
  -w "${TAXA_DIR}/${SAMPLE}_bracken_species_report.txt" \
  -r S \
  -l 150 \
  -t 10

echo "[STEP 5c] Re-estimating genus-level abundance with Bracken..."
bracken \
  -d "${KRAKEN2_DB}" \
  -i "${TAXA_DIR}/${SAMPLE}_kraken2_report.txt" \
  -o "${TAXA_DIR}/${SAMPLE}_bracken_genus.txt" \
  -w "${TAXA_DIR}/${SAMPLE}_bracken_genus_report.txt" \
  -r G \
  -l 150 \
  -t 10

echo "[STEP 5] Taxonomic profiling complete. Reports in ${TAXA_DIR}"

# ==============================================================================
# STEP 6: HUMAnN 3 Functional Profiling
# ==============================================================================
echo "[STEP 6] Running HUMAnN 3 for functional gene & pathway profiling..."
# Concatenate non-human paired reads for HUMAnN3 input
cat "${DECONTAM_DIR}/${SAMPLE}_nonhuman_1.fastq.gz" \
    "${DECONTAM_DIR}/${SAMPLE}_nonhuman_2.fastq.gz" \
    > "${DECONTAM_DIR}/${SAMPLE}_nonhuman_merged.fastq.gz"

humann \
  --input "${DECONTAM_DIR}/${SAMPLE}_nonhuman_merged.fastq.gz" \
  --output "${FUNC_DIR}" \
  --threads "${THREADS}" \
  --nucleotide-database "${HUMANN_NUC_DB}" \
  --protein-database "${HUMANN_PROT_DB}"

# Normalize gene families and pathway abundances to relative abundance (CPM: Copies Per Million)
humann_renorm_table \
  --input "${FUNC_DIR}/${SAMPLE}_nonhuman_merged_genefamilies.tsv" \
  --output "${FUNC_DIR}/${SAMPLE}_genefamilies_cpm.tsv" \
  --units cpm

humann_renorm_table \
  --input "${FUNC_DIR}/${SAMPLE}_nonhuman_merged_pathabundance.tsv" \
  --output "${FUNC_DIR}/${SAMPLE}_pathabundance_cpm.tsv" \
  --units cpm

echo "[STEP 6] Functional profiling complete. CPM-normalized tables in ${FUNC_DIR}"

# ==============================================================================
# STEP 7: MaAsLin 3 Cohort Multivariate Association Analysis
# ==============================================================================
echo "[STEP 7] Multi-sample multivariate association analysis using MaAsLin 3..."
echo "Note: MaAsLin3 requires cohort metadata across N samples to test for statistical"
echo "associations (e.g. CRC tumor vs control, age, sex). Running R script..."

Rscript run_maaslin3_analysis.R

echo "======================================================================"
echo "Shotgun Metagenomics Pipeline Execution Complete!"
echo "======================================================================"
