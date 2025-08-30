# Direct APOE Analysis - Single Cell RNA-seq
# Adapted from 1000 Genomes APOE workflow for RNA-seq data
# Focus: APOE gene expression patterns across IDH genotypes

library(dplyr)
library(ggplot2)

# Create output directories
dir.create("results/APOE_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/APOE_analysis", recursive = TRUE, showWarnings = FALSE)

cat("=== DIRECT APOE ANALYSIS - SINGLE CELL RNA-SEQ ===\n")

# Load metadata
metadata <- read.csv("data/metadata/GSE222520_phenotype_data.csv", stringsAsFactors = FALSE)

# Create sample info with proper column names
sample_info <- data.frame(
  sample_id = gsub(", replicate1, scRNAseq", "", metadata$title),
  geo_accession = metadata$geo_accession,
  tissue = metadata$tissue.ch1,
  diagnosis = metadata$diagnosis.ch1,
  genotype = metadata$genotype.ch1,
  primary_recurrent = metadata$primary.recurrent.ch1,
  stringsAsFactors = FALSE
)

# Define groups (similar to population stratification in 1000 Genomes)
sample_info$group <- case_when(
  grepl("^NGB", sample_info$sample_id) ~ "Normal",
  grepl("^IMP", sample_info$sample_id) ~ "Primary_IDH_Mutant",
  grepl("^IWP", sample_info$sample_id) ~ "Primary_IDH_WildType",
  grepl("^IMR", sample_info$sample_id) ~ "Recurrent_IDH_Mutant",
  grepl("^IWR", sample_info$sample_id) ~ "Recurrent_IDH_WildType"
)

sample_info$apoe_group <- case_when(
  sample_info$genotype == "IDH Mutant" ~ "IDH_Mutant",
  sample_info$genotype == "IDH Wild Type" ~ "IDH_WildType",
  TRUE ~ "Normal"
)

# APOE gene coordinates (from your 1000 Genomes analysis)
# APOE: chr19:45411941-45412079 (rs429358, rs7412 region)
# In RNA-seq, we look for APOE gene expression

# Function to extract APOE expression directly
extract_apoe_expression <- function(sample_dir) {
  sample_name <- gsub("GSM[0-9]+_", "", sample_dir)
  
  # Find matrix directory
  matrix_path <- file.path("data/raw/extracted", sample_dir)
  subdirs <- list.dirs(matrix_path, recursive = FALSE)
  if (length(subdirs) > 0) {
    matrix_path <- file.path(subdirs[1], "filtered_feature_bc_matrix")
  }
  
  if (!dir.exists(matrix_path)) {
    return(NULL)
  }
  
  # Read features file
  features_file <- file.path(matrix_path, "features.tsv.gz")
  if (!file.exists(features_file)) {
    return(NULL)
  }
  
  tryCatch({
    features <- readLines(gzfile(features_file))
    
    # Find APOE gene index
    apoe_index <- NULL
    for (i in 1:length(features)) {
      parts <- strsplit(features[i], "\t")[[1]]
      if (parts[2] == "APOE") {
        apoe_index <- i
        break
      }
    }
    
    if (!is.null(apoe_index)) {
      # Read matrix file
      matrix_file <- file.path(matrix_path, "matrix.mtx.gz")
      matrix_lines <- readLines(gzfile(matrix_file))
      
      # Parse header
      header <- strsplit(matrix_lines[1], " ")[[1]]
      n_genes <- as.numeric(header[1])
      n_cells <- as.numeric(header[2])
      
      # Parse matrix entries
      entries <- matrix_lines[2:length(matrix_lines)]
      parsed_entries <- do.call(rbind, lapply(entries, function(x) {
        parts <- strsplit(x, " ")[[1]]
        c(as.numeric(parts[1]), as.numeric(parts[2]), as.numeric(parts[3]))
      }))
      
      # Get APOE expression entries
      apoe_entries <- parsed_entries[parsed_entries[, 1] == apoe_index, , drop = FALSE]
      
      if (nrow(apoe_entries) > 0) {
        # Calculate APOE expression metrics
        mean_expr <- mean(apoe_entries[, 3])
        median_expr <- median(apoe_entries[, 3])
        max_expr <- max(apoe_entries[, 3])
        n_expressing_cells <- nrow(apoe_entries)
        pct_expressing <- n_expressing_cells / n_cells * 100
        
        # Calculate expression percentiles (similar to allele frequencies)
        expr_percentiles <- quantile(apoe_entries[, 3], probs = c(0.25, 0.5, 0.75, 0.9, 0.95))
        
        return(data.frame(
          sample = sample_name,
          mean_expression = mean_expr,
          median_expression = median_expr,
          max_expression = max_expr,
          n_expressing_cells = n_expressing_cells,
          total_cells = n_cells,
          pct_expressing = pct_expressing,
          expr_q25 = expr_percentiles[1],
          expr_q50 = expr_percentiles[2],
          expr_q75 = expr_percentiles[3],
          expr_q90 = expr_percentiles[4],
          expr_q95 = expr_percentiles[5],
          stringsAsFactors = FALSE
        ))
      } else {
        return(data.frame(
          sample = sample_name,
          mean_expression = 0,
          median_expression = 0,
          max_expression = 0,
          n_expressing_cells = 0,
          total_cells = n_cells,
          pct_expressing = 0,
          expr_q25 = 0,
          expr_q50 = 0,
          expr_q75 = 0,
          expr_q90 = 0,
          expr_q95 = 0,
          stringsAsFactors = FALSE
        ))
      }
    } else {
      return(data.frame(
        sample = sample_name,
        mean_expression = 0,
        median_expression = 0,
        max_expression = 0,
        n_expressing_cells = 0,
        total_cells = 0,
        pct_expressing = 0,
        expr_q25 = 0,
        expr_q50 = 0,
        expr_q75 = 0,
        expr_q90 = 0,
        expr_q95 = 0,
        stringsAsFactors = FALSE
      ))
    }
  }, error = function(e) {
    return(data.frame(
      sample = sample_name,
      mean_expression = 0,
      median_expression = 0,
      max_expression = 0,
      n_expressing_cells = 0,
      total_cells = 0,
      pct_expressing = 0,
      expr_q25 = 0,
      expr_q50 = 0,
      expr_q75 = 0,
      expr_q90 = 0,
      expr_q95 = 0,
      stringsAsFactors = FALSE
    ))
  })
}

# Get sample directories
sample_dirs <- list.dirs("data/raw/extracted", recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[!grepl("^\\.", sample_dirs)]

cat("Analyzing APOE expression in", length(sample_dirs), "samples...\n")

# Extract APOE data from all samples
all_results <- lapply(sample_dirs, extract_apoe_expression)
all_results <- do.call(rbind, all_results)

# Combine with metadata
if (!is.null(all_results) && nrow(all_results) > 0) {
  all_results <- merge(all_results, sample_info, by.x = "sample", by.y = "sample_id", all.x = TRUE)
  
  # Save results
  write.csv(all_results, "results/APOE_analysis/direct_APOE_results.csv", row.names = FALSE)
  
  # Create summary (similar to 1000 Genomes frequency analysis)
  cat("\n=== APOE EXPRESSION SUMMARY ===\n")
  cat("Total samples analyzed:", length(unique(all_results$sample)), "\n")
  
  # Overall APOE expression summary
  overall_summary <- all_results %>%
    summarise(
      n_samples = n(),
      mean_expression = mean(mean_expression, na.rm = TRUE),
      median_expression = median(median_expression, na.rm = TRUE),
      mean_pct_expressing = mean(pct_expressing, na.rm = TRUE),
      max_expression = max(max_expression, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\nOverall APOE Expression Summary:\n")
  print(overall_summary)
  
  # Group-wise summary (similar to population stratification)
  group_summary <- all_results %>%
    group_by(apoe_group) %>%
    summarise(
      n_samples = n(),
      mean_expression = mean(mean_expression, na.rm = TRUE),
      median_expression = median(median_expression, na.rm = TRUE),
      mean_pct_expressing = mean(pct_expressing, na.rm = TRUE),
      max_expression = max(max_expression, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\nAPOE Expression by IDH Genotype Group:\n")
  print(group_summary)
  
  # Diagnosis-wise summary
  diagnosis_summary <- all_results %>%
    group_by(diagnosis) %>%
    summarise(
      n_samples = n(),
      mean_expression = mean(mean_expression, na.rm = TRUE),
      mean_pct_expressing = mean(pct_expressing, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\nAPOE Expression by Diagnosis:\n")
  print(diagnosis_summary)
  
  # Create visualizations (similar to 1000 Genomes plots)
  if (nrow(all_results) > 0) {
    # APOE expression distribution by group
    p1 <- ggplot(all_results, aes(x = apoe_group, y = mean_expression, fill = apoe_group)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.6) +
      theme_minimal() +
      labs(title = "APOE Expression by IDH Genotype Group",
           x = "IDH Genotype Group", y = "Mean Expression", fill = "Group")
    
    ggsave("figures/APOE_analysis/direct_APOE_expression_by_group.pdf", p1, width = 10, height = 6)
    
    # APOE expression percentage by group
    p2 <- ggplot(all_results, aes(x = apoe_group, y = pct_expressing, fill = apoe_group)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.6) +
      theme_minimal() +
      labs(title = "APOE Expression Percentage by IDH Genotype Group",
           x = "IDH Genotype Group", y = "Percentage of Expressing Cells", fill = "Group")
    
    ggsave("figures/APOE_analysis/direct_APOE_percentage_by_group.pdf", p2, width = 10, height = 6)
    
    # APOE expression by diagnosis
    p3 <- ggplot(all_results, aes(x = diagnosis, y = mean_expression, fill = apoe_group)) +
      geom_boxplot() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "APOE Expression by Diagnosis and IDH Status",
           x = "Diagnosis", y = "Mean Expression", fill = "IDH Group")
    
    ggsave("figures/APOE_analysis/direct_APOE_by_diagnosis.pdf", p3, width = 12, height = 6)
    
    # Expression distribution histogram
    p4 <- ggplot(all_results, aes(x = mean_expression, fill = apoe_group)) +
      geom_histogram(bins = 20, alpha = 0.7) +
      facet_wrap(~apoe_group) +
      theme_minimal() +
      labs(title = "APOE Expression Distribution by IDH Group",
           x = "Mean Expression", y = "Count", fill = "IDH Group")
    
    ggsave("figures/APOE_analysis/direct_APOE_distribution.pdf", p4, width = 12, height = 8)
  }
  
  cat("\nResults saved to: results/APOE_analysis/direct_APOE_results.csv\n")
  cat("Figures saved to: figures/APOE_analysis/\n")
  
} else {
  cat("No results obtained from APOE analysis.\n")
}

cat("\n=== DIRECT APOE ANALYSIS COMPLETE ===\n")
cat("Note: This analysis focuses on APOE gene expression patterns.\n")
cat("For APOE isoforms (ε2/ε3/ε4), germline SNP data (rs429358, rs7412) is required.\n")
