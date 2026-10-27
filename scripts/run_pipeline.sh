#!/usr/bin/env bash
# ==============================================================================
# Master Pipeline: run_pipeline.sh
# End-to-End Shotgun Metagenomics Workflow (Healthy vs CRC)
# FastQC -> fastp -> Bowtie2 -> Kraken2/Bracken -> HUMAnN3 -> MaAsLin3
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "${SCRIPT_DIR}")"
cd "${BASE_DIR}"

echo "========================================================================"
echo " Starting Shotgun Metagenomics Pipeline: CRC vs Healthy Cohort"
echo " Accessions: SRX20110651, SRX20110652 (CRC) | SRX20110653, SRX20110654 (Healthy)"
echo " Start Time: $(date)"
echo "========================================================================"

# Step 0: Download Raw Data from NCBI SRA
echo -e "\n>>> STEP 0: Downloading SRA Datasets..."
bash "${SCRIPT_DIR}/00_download_sra.sh"

# Step 1: Quality Control (Raw FastQC)
echo -e "\n>>> STEP 1: FastQC Quality Control..."
bash "${SCRIPT_DIR}/01_fastqc_raw.sh"

# Step 2: Quality Filtering & Adapter Trimming (fastp)
echo -e "\n>>> STEP 2: fastp Filtering and Trimming..."
bash "${SCRIPT_DIR}/02_fastp_trim.sh"

# Step 3: Host Decontamination (Bowtie2 - Human GRCh38 Removal)
echo -e "\n>>> STEP 3: Bowtie2 Host Decontamination..."
bash "${SCRIPT_DIR}/03_bowtie2_decontam.sh"

# Step 4: Taxonomic Profiling (Kraken2 & Bracken)
echo -e "\n>>> STEP 4: Kraken2 Taxonomic Profiling & Bracken Re-estimation..."
bash "${SCRIPT_DIR}/04_kraken2_bracken.sh"

# Step 5: Functional Profiling (HUMAnN3)
echo -e "\n>>> STEP 5: HUMAnN3 Functional & Metabolic Pathway Profiling..."
bash "${SCRIPT_DIR}/05_humann3_profiling.sh"

# Step 6: Multivariable Association & Statistical Testing (MaAsLin3)
echo -e "\n>>> STEP 6: MaAsLin3 Differential Abundance Analysis..."
Rscript "${SCRIPT_DIR}/06_maaslin3_analysis.R" \
    --taxa "04_bracken_output/merged_bracken_species_abundance.txt" \
    --pathways "05_humann3_joined/cohort_pathabundance_cpm_unstratified.tsv" \
    --metadata "metadata.tsv" \
    --outdir "06_maaslin3_results"

echo -e "\n========================================================================"
echo " Shotgun Metagenomics Pipeline Finished Successfully!"
echo " End Time: $(date)"
echo "========================================================================"
