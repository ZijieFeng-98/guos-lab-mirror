# =============================================================================
# COMPREHENSIVE APOE ISOFORM ANALYSIS WORKFLOW
# =============================================================================
# 
# This script implements a complete APOE isoform analysis pipeline that can be
# adapted for different genomic datasets including 1000 Genomes and TCGA.
#
# APOE LOGIC OVERVIEW:
# ===================
# APOE (Apolipoprotein E) has three major isoforms: ε2, ε3, and ε4, determined
# by two SNPs on chromosome 19:
# 
# 1. rs429358 (chr19:45411941, GRCh37/hg19)
#    - Reference allele: T
#    - Alternate allele: C
# 
# 2. rs7412 (chr19:45412079, GRCh37/hg19)  
#    - Reference allele: C
#    - Alternate allele: T
#
# HAPLOTYPE ASSIGNMENT LOGIC:
# ===========================
# From phased genotypes (e.g., 0|1, 1|0), we construct haplotypes:
# - Haplotype 1: rs429358_allele1 + rs7412_allele1
# - Haplotype 2: rs429358_allele2 + rs7412_allele2
#
# APOE isoform assignment based on haplotype:
# - ε2: rs429358=T, rs7412=T (T-T haplotype)
# - ε3: rs429358=T, rs7412=C (T-C haplotype) 
# - ε4: rs429358=C, rs7412=C (C-C haplotype)
# - ε1: rs429358=C, rs7412=T (C-T haplotype, non-canonical)
#
# GENOTYPE ASSIGNMENT:
# ====================
# Individual genotypes are the combination of two haplotypes:
# - ε2/ε2, ε2/ε3, ε2/ε4, ε3/ε3, ε3/ε4, ε4/ε4, etc.
# - Sorted alphabetically (e.g., ε2/ε4 not ε4/ε2)
#
# CARRIER STATUS:
# ===============
# - ε4 carrier: Has at least one ε4 allele
# - ε2 carrier: Has at least one ε2 allele
#
# CLINICAL RELEVANCE:
# ===================
# - ε4: Associated with increased Alzheimer's disease risk, cardiovascular disease
# - ε2: Associated with decreased Alzheimer's risk, but increased cardiovascular risk
# - ε3: Most common, considered "neutral"
# - Population differences: ε4 higher in African populations, ε2 higher in Asian populations
#
# ADAPTATION FOR TCGA:
# ===================
# For TCGA cancer data analysis:
# 1. Use the same SNP coordinates (rs429358, rs7412)
# 2. Apply the same haplotype assignment logic
# 3. Stratify by cancer type instead of population
# 4. Compare APOE frequencies between cancer types and controls
# 5. Test associations between APOE isoforms and cancer outcomes
#
# =============================================================================

# Load required libraries
library(vcfR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(purrr)
library(forcats)

# =============================================================================
# FUNCTION DEFINITIONS
# =============================================================================

#' Verify SNP coordinates in GRCh37/hg19
#' @return Data frame with SNP information
verify_snp_coordinates <- function() {
  cat("=== Verifying APOE SNP Coordinates (GRCh37/hg19) ===\n")
  
  snp_info <- data.frame(
    rsid = c("rs429358", "rs7412"),
    chr = c(19, 19),
    pos = c(45411941, 45412079),  # GRCh37 coordinates
    ref = c("T", "C"),
    alt = c("C", "T"),
    stringsAsFactors = FALSE
  )
  
  cat("APOE SNPs:\n")
  print(snp_info)
  cat("\n")
  
  return(snp_info)
}

#' Download 1000 Genomes data if not present
#' @param force_download Logical, force re-download if TRUE
download_1000g_data <- function(force_download = FALSE) {
  cat("=== Downloading 1000 Genomes Data ===\n")
  
  # Create data directory
  if (!dir.exists("data")) dir.create("data")
  
  # URLs for 1000 Genomes Phase 3 data
  vcf_url <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  sample_url <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
  
  # File paths
  vcf_file <- "data/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  sample_file <- "data/integrated_call_samples_v3.20130502.ALL.panel"
  
  # Download VCF if not present or forced
  if (!file.exists(vcf_file) || force_download) {
    cat("Downloading chromosome 19 VCF file...\n")
    download.file(vcf_url, vcf_file, method = "auto")
    
    # Download tabix index
    tabix_url <- paste0(vcf_url, ".tbi")
    tabix_file <- paste0(vcf_file, ".tbi")
    download.file(tabix_url, tabix_file, method = "auto")
  } else {
    cat("VCF file already exists\n")
  }
  
  # Download sample metadata if not present or forced
  if (!file.exists(sample_file) || force_download) {
    cat("Downloading sample metadata...\n")
    download.file(sample_url, sample_file, method = "auto")
  } else {
    cat("Sample metadata already exists\n")
  }
  
  cat("Data download complete\n\n")
}

#' Extract specific SNPs from VCF file
#' @param vcf_file Path to VCF file
#' @param snp_info SNP information data frame
#' @return Path to extracted VCF file
extract_snps_from_vcf <- function(vcf_file, snp_info) {
  cat("=== Extracting APOE SNPs from VCF ===\n")
  
  # Create temporary file for extraction
  temp_vcf <- tempfile(fileext = ".vcf")
  
  # Extract by position since SNPs may not have rsIDs in VCF
  pos_pattern <- paste(snp_info$pos, collapse = "|")
  grep_cmd <- paste("gunzip -c", vcf_file, "| grep -E '^#|", pos_pattern, "' >", temp_vcf)
  
  cat("Running extraction command...\n")
  system(grep_cmd)
  
  # Check if extraction worked
  if (!file.exists(temp_vcf) || file.size(temp_vcf) == 0) {
    stop("Failed to extract SNPs from VCF file")
  }
  
  cat("Extraction successful. File size:", file.size(temp_vcf), "bytes\n\n")
  return(temp_vcf)
}

#' Load and process VCF data
#' @param vcf_file Path to VCF file
#' @return Processed VCF object
load_vcf_data <- function(vcf_file) {
  cat("=== Loading VCF Data ===\n")
  
  # Read VCF file
  vcf <- read.vcfR(vcf_file)
  
  # Check what we have
  cat("VCF info:\n")
  cat("  Variants:", nrow(vcf), "\n")
  cat("  Samples:", ncol(vcf), "\n")
  
  if (nrow(vcf) == 0) {
    stop("No variants found in VCF file")
  }
  
  # Extract SNP information
  snp_data <- getFIX(vcf)
  cat("SNP positions:", paste(snp_data[, "POS"], collapse = ", "), "\n\n")
  
  return(vcf)
}

#' Extract genotypes and perform quality control
#' @param vcf VCF object
#' @return Genotype data frame
extract_genotypes <- function(vcf) {
  cat("=== Extracting Genotypes ===\n")
  
  # Extract genotype matrix
  gt_matrix <- extract.gt(vcf)
  cat("Genotype matrix dimensions:", dim(gt_matrix), "\n")
  
  # Get sample names
  sample_names <- colnames(gt_matrix)
  
  # Create data frame with genotypes
  geno_df <- data.frame(
    sample_id = sample_names,
    rs429358 = gt_matrix[1, ],
    rs7412 = gt_matrix[2, ],
    stringsAsFactors = FALSE
  )
  
  # Basic QC: remove samples with missing genotypes
  initial_samples <- nrow(geno_df)
  geno_df <- geno_df[!is.na(geno_df$rs429358) & !is.na(geno_df$rs7412), ]
  final_samples <- nrow(geno_df)
  
  cat("QC: Removed", initial_samples - final_samples, "samples with missing genotypes\n")
  cat("Final sample count:", final_samples, "\n")
  
  # Verify phasing (should be 0|1, 1|0, etc.)
  phased_check <- all(grepl("\\|", geno_df$rs429358)) & all(grepl("\\|", geno_df$rs7412))
  cat("Phasing check:", ifelse(phased_check, "PASS", "FAIL"), "\n\n")
  
  return(geno_df)
}

#' Build haplotypes and assign APOE isoforms
#' @param geno_df Genotype data frame
#' @return Haplotype data frame
build_haplotypes <- function(geno_df) {
  cat("=== Building Haplotypes and Assigning APOE Isoforms ===\n")
  
  # Function to split phased genotypes
  split_phased <- function(gt_string) {
    if (is.na(gt_string) || gt_string == "") return(c(NA, NA))
    alleles <- str_split(gt_string, "\\|")[[1]]
    return(as.numeric(alleles))
  }
  
  # Split genotypes for both SNPs
  rs429358_split <- t(sapply(geno_df$rs429358, split_phased))
  rs7412_split <- t(sapply(geno_df$rs7412, split_phased))
  
  # Create haplotype data frame
  haplo_df <- data.frame(
    sample_id = geno_df$sample_id,
    hap1_rs429358 = rs429358_split[, 1],
    hap1_rs7412 = rs7412_split[, 1],
    hap2_rs429358 = rs429358_split[, 2],
    hap2_rs7412 = rs7412_split[, 2],
    stringsAsFactors = FALSE
  )
  
  # Function to assign APOE haplotype
  assign_apoe_haplotype <- function(rs429358_allele, rs7412_allele) {
    # Convert numeric alleles to actual bases (0=REF, 1=ALT)
    # For rs429358: REF=T, ALT=C
    # For rs7412: REF=C, ALT=T
    
    rs429358_base <- ifelse(rs429358_allele == 0, "T", "C")
    rs7412_base <- ifelse(rs7412_allele == 0, "C", "T")
    
    # Assign APOE haplotype based on logic
    if (rs429358_base == "T" && rs7412_base == "T") return("ε2")
    if (rs429358_base == "T" && rs7412_base == "C") return("ε3")
    if (rs429358_base == "C" && rs7412_base == "C") return("ε4")
    if (rs429358_base == "C" && rs7412_base == "T") return("ε1") # Non-canonical
    
    return("unknown")
  }
  
  # Assign haplotypes
  haplo_df$hap1_apoe <- mapply(assign_apoe_haplotype, 
                               haplo_df$hap1_rs429358, 
                               haplo_df$hap1_rs7412)
  haplo_df$hap2_apoe <- mapply(assign_apoe_haplotype, 
                               haplo_df$hap2_rs429358, 
                               haplo_df$hap2_rs7412)
  
  # Create genotype (sorted haplotypes)
  haplo_df$genotype <- paste(
    pmin(haplo_df$hap1_apoe, haplo_df$hap2_apoe),
    pmax(haplo_df$hap1_apoe, haplo_df$hap2_apoe),
    sep = "/"
  )
  
  # Add carrier status
  haplo_df$carrier_e4 <- haplo_df$hap1_apoe == "ε4" | haplo_df$hap2_apoe == "ε4"
  haplo_df$carrier_e2 <- haplo_df$hap1_apoe == "ε2" | haplo_df$hap2_apoe == "ε2"
  
  cat("Haplotype assignment completed\n")
  cat("Sample genotypes (first 10):\n")
  print(head(haplo_df[, c("sample_id", "hap1_apoe", "hap2_apoe", "genotype")], 10))
  cat("\n")
  
  return(haplo_df)
}

#' Load sample metadata
#' @param sample_file Path to sample metadata file
#' @return Metadata data frame
load_sample_metadata <- function(sample_file) {
  cat("=== Loading Sample Metadata ===\n")
  
  metadata <- read_tsv(sample_file, col_names = c("sample_id", "pop", "superpop", "sex"))
  cat("Loaded metadata for", nrow(metadata), "samples\n\n")
  
  return(metadata)
}

#' Merge genotype and metadata
#' @param haplo_df Haplotype data frame
#' @param metadata Metadata data frame
#' @return Merged data frame
merge_data <- function(haplo_df, metadata) {
  cat("=== Merging Genotype and Metadata ===\n")
  
  merged_df <- left_join(haplo_df, metadata, by = "sample_id")
  cat("Merged data dimensions:", dim(merged_df), "\n\n")
  
  return(merged_df)
}

#' Compute summary statistics
#' @param merged_df Merged data frame
#' @return List of summary statistics
compute_summaries <- function(merged_df) {
  cat("=== Computing Summary Statistics ===\n")
  
  # Overall haplotype frequencies
  hap_freq <- merged_df %>%
    select(hap1_apoe, hap2_apoe) %>%
    pivot_longer(cols = everything(), names_to = "haplotype", values_to = "allele") %>%
    count(allele) %>%
    mutate(frequency = n / sum(n)) %>%
    arrange(desc(frequency))
  
  # Genotype frequencies
  geno_freq <- merged_df %>%
    count(genotype) %>%
    mutate(frequency = n / sum(n)) %>%
    arrange(desc(frequency))
  
  # Carrier rates
  carrier_e4_rate <- mean(merged_df$carrier_e4, na.rm = TRUE)
  carrier_e2_rate <- mean(merged_df$carrier_e2, na.rm = TRUE)
  
  # By superpopulation
  superpop_summary <- merged_df %>%
    group_by(superpop) %>%
    summarise(
      n_samples = n(),
      e4_carrier_rate = mean(carrier_e4, na.rm = TRUE),
      e2_carrier_rate = mean(carrier_e2, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Print results
  cat("Overall haplotype frequencies:\n")
  print(hap_freq)
  cat("\nGenotype frequencies:\n")
  print(geno_freq)
  cat("\nε4 carrier rate:", round(carrier_e4_rate * 100, 2), "%\n")
  cat("ε2 carrier rate:", round(carrier_e2_rate * 100, 2), "%\n")
  cat("\nSummary by superpopulation:\n")
  print(superpop_summary)
  cat("\n")
  
  return(list(
    haplotype_frequencies = hap_freq,
    genotype_frequencies = geno_freq,
    carrier_e4_rate = carrier_e4_rate,
    carrier_e2_rate = carrier_e2_rate,
    superpopulation_summary = superpop_summary
  ))
}

#' Perform statistical tests
#' @param merged_df Merged data frame
#' @return Statistical test results
perform_statistical_tests <- function(merged_df) {
  cat("=== Performing Statistical Tests ===\n")
  
  # Chi-square test for ε4 carrier vs non-carrier by superpopulation
  e4_table <- table(merged_df$superpop, merged_df$carrier_e4)
  e4_chi <- chisq.test(e4_table)
  
  # Chi-square test for ε2 carrier vs non-carrier by superpopulation
  e2_table <- table(merged_df$superpop, merged_df$carrier_e2)
  e2_chi <- chisq.test(e2_table)
  
  # Calculate odds ratios for ε4 carrier between population pairs
  e4_rates <- merged_df %>%
    group_by(superpop) %>%
    summarise(e4_rate = mean(carrier_e4, na.rm = TRUE), .groups = "drop")
  
  cat("Chi-square test for ε4 carrier by superpopulation:\n")
  print(e4_chi)
  cat("\nChi-square test for ε2 carrier by superpopulation:\n")
  print(e2_chi)
  cat("\nε4 carrier rates by population:\n")
  print(e4_rates)
  cat("\n")
  
  return(list(
    e4_chi_square = e4_chi,
    e2_chi_square = e2_chi,
    e4_rates_by_pop = e4_rates
  ))
}

#' Create visualizations
#' @param summaries Summary statistics
#' @param stats_results Statistical test results
create_visualizations <- function(summaries, stats_results) {
  cat("=== Creating Visualizations ===\n")
  
  # Create output directory
  if (!dir.exists("output")) dir.create("output")
  
  # 1. Haplotype frequency bar plot
  p1 <- ggplot(summaries$haplotype_frequencies, aes(x = allele, y = frequency)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = sprintf("%.1f%%", frequency * 100)), vjust = -0.5) +
    labs(title = "APOE Haplotype Frequencies", x = "Haplotype", y = "Frequency") +
    theme_minimal()
  
  ggsave("output/haplotype_frequencies.png", p1, width = 8, height = 6)
  
  # 2. Genotype frequency bar plot
  p2 <- ggplot(summaries$genotype_frequencies, aes(x = genotype, y = frequency)) +
    geom_bar(stat = "identity", fill = "darkgreen") +
    geom_text(aes(label = sprintf("%.1f%%", frequency * 100)), vjust = -0.5, angle = 45) +
    labs(title = "APOE Genotype Frequencies", x = "Genotype", y = "Frequency") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("output/genotype_frequencies.png", p2, width = 10, height = 6)
  
  # 3. ε4 carrier rates by population
  p3 <- ggplot(summaries$superpopulation_summary, aes(x = superpop, y = e4_carrier_rate)) +
    geom_bar(stat = "identity", fill = "red") +
    geom_text(aes(label = sprintf("%.1f%%", e4_carrier_rate * 100)), vjust = -0.5) +
    labs(title = "ε4 Carrier Rates by Population", x = "Population", y = "ε4 Carrier Rate") +
    theme_minimal()
  
  ggsave("output/e4_carrier_by_population.png", p3, width = 8, height = 6)
  
  cat("Visualizations saved to output/ directory\n\n")
}

#' Export results
#' @param merged_df Merged data frame
#' @param summaries Summary statistics
#' @param stats_results Statistical test results
export_results <- function(merged_df, summaries, stats_results) {
  cat("=== Exporting Results ===\n")
  
  # Create output directory
  if (!dir.exists("output")) dir.create("output")
  
  # Save data files
  write_csv(merged_df, "output/apoe_genotypes.csv")
  write_csv(summaries$haplotype_frequencies, "output/haplotype_frequencies.csv")
  write_csv(summaries$genotype_frequencies, "output/genotype_frequencies.csv")
  write_csv(summaries$superpopulation_summary, "output/superpopulation_summary.csv")
  
  # Save statistical test results
  capture.output(summaries, file = "output/statistical_tests.txt")
  
  cat("Results exported to output/ directory\n\n")
}

#' Main analysis function
#' @param force_download Logical, force re-download if TRUE
#' @return List of all results
run_apoe_analysis <- function(force_download = FALSE) {
  cat("=== COMPREHENSIVE APOE ANALYSIS WORKFLOW ===\n\n")
  
  # Step 1: Verify SNP coordinates
  snp_info <- verify_snp_coordinates()
  
  # Step 2: Download data
  download_1000g_data(force_download)
  
  # Step 3: Extract SNPs
  vcf_file <- "data/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  extracted_vcf <- extract_snps_from_vcf(vcf_file, snp_info)
  
  # Step 4: Load VCF data
  vcf <- load_vcf_data(extracted_vcf)
  
  # Step 5: Extract genotypes
  geno_df <- extract_genotypes(vcf)
  
  # Step 6: Build haplotypes
  haplo_df <- build_haplotypes(geno_df)
  
  # Step 7: Load metadata
  sample_file <- "data/integrated_call_samples_v3.20130502.ALL.panel"
  metadata <- load_sample_metadata(sample_file)
  
  # Step 8: Merge data
  merged_df <- merge_data(haplo_df, metadata)
  
  # Step 9: Compute summaries
  summaries <- compute_summaries(merged_df)
  
  # Step 10: Perform statistical tests
  stats_results <- perform_statistical_tests(merged_df)
  
  # Step 11: Create visualizations
  create_visualizations(summaries, stats_results)
  
  # Step 12: Export results
  export_results(merged_df, summaries, stats_results)
  
  # Clean up temporary file
  unlink(extracted_vcf)
  
  cat("=== ANALYSIS COMPLETE ===\n")
  cat("All results saved to output/ directory\n")
  
  return(list(
    merged_data = merged_df,
    summaries = summaries,
    statistical_tests = stats_results
  ))
}

# =============================================================================
# RUN ANALYSIS
# =============================================================================

# Run the complete analysis
results <- run_apoe_analysis(force_download = FALSE)

