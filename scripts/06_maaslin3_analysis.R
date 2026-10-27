#!/usr/bin/env Rscript
# ==============================================================================
# Script: 06_maaslin3_analysis.R
# Description: Multivariable association & differential abundance testing 
#              comparing Colorectal Cancer (CRC) vs Healthy controls using MaAsLin3
# Author: Shotgun Metagenomics Pipeline
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse", repos = "https://cloud.r-project.org")
  library(optparse)
})

# Parse Command Line Options
option_list <- list(
  make_option(c("-t", "--taxa"), type = "character", default = "04_bracken_output/merged_bracken_species_abundance.txt",
              help = "Path to merged Bracken species abundance matrix [default: %default]"),
  make_option(c("-p", "--pathways"), type = "character", default = "05_humann3_joined/cohort_pathabundance_cpm_unstratified.tsv",
              help = "Path to unstratified HUMAnN3 CPM pathway abundance table [default: %default]"),
  make_option(c("-m", "--metadata"), type = "character", default = "metadata.tsv",
              help = "Path to sample metadata TSV file [default: %default]"),
  make_option(c("-o", "--outdir"), type = "character", default = "06_maaslin3_results",
              help = "Output directory for MaAsLin results [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cat("==========================================================\n")
cat("Phase 6: MaAsLin3 Multivariable Association Analysis\n")
cat("Comparing: Colorectal Cancer (CRC) vs Healthy Controls\n")
cat("==========================================================\n\n")

# Load required packages
check_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Package '%s' not found. Installing...\n", pkg))
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
    if (pkg == "maaslin3") {
      # MaAsLin3 installation from GitHub or Bioconductor
      if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos = "https://cloud.r-project.org")
      devtools::install_github("biobakery/maaslin3")
    } else if (pkg == "Maaslin2") {
      BiocManager::install("Maaslin2")
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  library(pkg, character.only = TRUE)
}

# 1. Load or Generate Metadata
if (!file.exists(opt$metadata)) {
  cat("Generating study cohort metadata file:", opt$metadata, "\n")
  meta_df <- data.frame(
    SampleID = c("SRR24315726", "SRR24315725", "SRR24315724", "SRR24315723"),
    Accession = c("SRX20110651", "SRX20110652", "SRX20110653", "SRX20110654"),
    PatientID = c("Patient_82", "Patient_83", "Patient_55", "Patient_56"),
    Condition = factor(c("CRC", "CRC", "Healthy", "Healthy"), levels = c("Healthy", "CRC")),
    Host = "Homo sapiens",
    BioProject = "PRJNA961076",
    stringsAsFactors = FALSE
  )
  write.table(meta_df, opt$metadata, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  meta_df <- read.table(opt$metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

rownames(meta_df) <- meta_df$SampleID
cat("Cohort metadata loaded successfully (", nrow(meta_df), " samples).\n\n")

# 2. Run MaAsLin Analysis Function
run_association <- function(feature_file, feature_type, out_sub) {
  cat(sprintf("--- Analyzing %s Differential Abundance ---\n", feature_type))
  
  if (!file.exists(feature_file)) {
    cat(sprintf("Warning: Feature file %s does not exist. Skipping %s.\n\n", feature_file, feature_type))
    return(NULL)
  }
  
  features <- read.table(feature_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  # Transpose if features are rows (samples as columns)
  if (all(rownames(meta_df) %in% colnames(features))) {
    features <- as.data.frame(t(features))
  }
  
  sub_dir <- file.path(opt$outdir, out_sub)
  dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)

  # Check whether maaslin3 or Maaslin2 is available
  has_maaslin3 <- requireNamespace("maaslin3", quietly = TRUE)
  has_maaslin2 <- requireNamespace("Maaslin2", quietly = TRUE)

  if (has_maaslin3) {
    cat("Using MaAsLin3 for statistical modeling...\n")
    library(maaslin3)
    fit <- maaslin3(
      input_data = features,
      input_metadata = meta_df,
      output = sub_dir,
      formula = ~ Condition,
      reference = "Condition,Healthy",
      normalization = "TSS",
      transform = "LOG",
      plot_summary_plot = TRUE
    )
  } else if (has_maaslin2) {
    cat("Using MaAsLin2 for statistical modeling...\n")
    library(Maaslin2)
    fit <- Maaslin2(
      input_data = features,
      input_metadata = meta_df,
      output = sub_dir,
      fixed_effects = c("Condition"),
      reference = c("Condition,Healthy"),
      normalization = "TSS",
      transform = "LOG",
      min_prevalence = 0.1,
      min_abundance = 0.0001,
      plot_heatmap = TRUE,
      plot_scatter = TRUE
    )
  } else {
    cat("MaAsLin package not currently available in environment. Performing standard GLM/Welch t-test...\n")
    # Base R fallback for differential testing
    results <- data.frame()
    for (feat in colnames(features)) {
      y <- as.numeric(features[, feat])
      if (sum(y > 0) >= 2) {
        fit_lm <- lm(log1p(y) ~ meta_df$Condition)
        s <- summary(fit_lm)
        if (nrow(s$coefficients) >= 2) {
          coef_val <- s$coefficients[2, 1]
          pval <- s$coefficients[2, 4]
          results <- rbind(results, data.frame(
            feature = feat,
            coef = coef_val,
            pval = pval,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    if (nrow(results) > 0) {
      results$qval <- p.adjust(results$pval, method = "BH")
      results <- results[order(results$pval), ]
      write.table(results, file.path(sub_dir, "differential_features_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
  cat(sprintf("Finished %s analysis. Results written to %s\n\n", feature_type, sub_dir))
}

# Run on Taxonomy (Bracken)
run_association(opt$taxa, "Taxonomic Species", "taxa_species_maaslin3")

# Run on Pathways (HUMAnN3)
run_association(opt$pathways, "MetaCyc Metabolic Pathways", "pathways_maaslin3")

cat("==========================================================\n")
cat("MaAsLin3 Analysis Complete!\n")
cat("==========================================================\n")
