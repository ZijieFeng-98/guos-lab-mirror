# APOE Best Practice Analysis - Single Cell RNA-seq
# Following best practices for APOE analysis without germline genotypes
# Focus: APOE expression patterns across conditions and cell types

library(dplyr)
library(ggplot2)
library(tidyr)

# Create output directories
dir.create("results/APOE_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/APOE_analysis", recursive = TRUE, showWarnings = FALSE)

cat("=== APOE BEST PRACTICE ANALYSIS ===\n")
cat("Data Type: Single-donor scRNA-seq samples\n")
cat("Limitation: No germline APOE genotype data (rs429358, rs7412)\n")
cat("Focus: APOE expression patterns across conditions\n\n")

# Load metadata
metadata <- read.csv("data/metadata/GSE222520_phenotype_data.csv", stringsAsFactors = FALSE)

# Create comprehensive sample info
sample_info <- data.frame(
  sample_id = gsub(", replicate1, scRNAseq", "", metadata$title),
  geo_accession = metadata$geo_accession,
  tissue = metadata$tissue.ch1,
  diagnosis = metadata$diagnosis.ch1,
  genotype = metadata$genotype.ch1,
  primary_recurrent = metadata$primary.recurrent.ch1,
  treatment = metadata$treatment.ch1,
  cell_type = metadata$cell.type.ch1,
  stringsAsFactors = FALSE
)

# Define comprehensive groups
sample_info$condition_group <- case_when(
  grepl("^NGB", sample_info$sample_id) ~ "Normal_Brain",
  grepl("^IMP", sample_info$sample_id) ~ "Primary_IDH_Mutant",
  grepl("^IWP", sample_info$sample_id) ~ "Primary_IDH_WildType", 
  grepl("^IMR", sample_info$sample_id) ~ "Recurrent_IDH_Mutant",
  grepl("^IWR", sample_info$sample_id) ~ "Recurrent_IDH_WildType"
)

sample_info$idh_status <- case_when(
  sample_info$genotype == "IDH Mutant" ~ "IDH_Mutant",
  sample_info$genotype == "IDH Wild Type" ~ "IDH_WildType",
  TRUE ~ "Normal"
)

sample_info$tumor_stage <- case_when(
  grepl("^NG", sample_info$sample_id) ~ "Normal",
  grepl("^IMP|^IWP", sample_info$sample_id) ~ "Primary",
  grepl("^IMR|^IWR", sample_info$sample_id) ~ "Recurrent"
)

# Function to extract APOE expression with comprehensive metrics
extract_apoe_comprehensive <- function(sample_dir) {
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
        # Comprehensive APOE expression metrics
        mean_expr <- mean(apoe_entries[, 3])
        median_expr <- median(apoe_entries[, 3])
        max_expr <- max(apoe_entries[, 3])
        n_expressing_cells <- nrow(apoe_entries)
        pct_expressing <- n_expressing_cells / n_cells * 100
        
        # Expression percentiles
        expr_percentiles <- quantile(apoe_entries[, 3], probs = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
        
        # Expression variance and distribution
        expr_variance <- var(apoe_entries[, 3])
        expr_cv <- sd(apoe_entries[, 3]) / mean_expr  # coefficient of variation
        
        # Expression categories (similar to genotype categories)
        high_expr_cells <- sum(apoe_entries[, 3] > quantile(apoe_entries[, 3], 0.75))
        low_expr_cells <- sum(apoe_entries[, 3] < quantile(apoe_entries[, 3], 0.25))
        
        return(data.frame(
          sample = sample_name,
          donor_id = sample_name,  # Each sample = one donor
          mean_expression = mean_expr,
          median_expression = median_expr,
          max_expression = max_expr,
          n_expressing_cells = n_expressing_cells,
          total_cells = n_cells,
          pct_expressing = pct_expressing,
          expr_variance = expr_variance,
          expr_cv = expr_cv,
          expr_q10 = expr_percentiles[1],
          expr_q25 = expr_percentiles[2],
          expr_q50 = expr_percentiles[3],
          expr_q75 = expr_percentiles[4],
          expr_q90 = expr_percentiles[5],
          expr_q95 = expr_percentiles[6],
          high_expr_cells = high_expr_cells,
          low_expr_cells = low_expr_cells,
          pct_high_expr = high_expr_cells / n_cells * 100,
          pct_low_expr = low_expr_cells / n_cells * 100,
          stringsAsFactors = FALSE
        ))
      } else {
        return(data.frame(
          sample = sample_name,
          donor_id = sample_name,
          mean_expression = 0,
          median_expression = 0,
          max_expression = 0,
          n_expressing_cells = 0,
          total_cells = n_cells,
          pct_expressing = 0,
          expr_variance = 0,
          expr_cv = 0,
          expr_q10 = 0,
          expr_q25 = 0,
          expr_q50 = 0,
          expr_q75 = 0,
          expr_q90 = 0,
          expr_q95 = 0,
          high_expr_cells = 0,
          low_expr_cells = 0,
          pct_high_expr = 0,
          pct_low_expr = 0,
          stringsAsFactors = FALSE
        ))
      }
    } else {
      return(data.frame(
        sample = sample_name,
        donor_id = sample_name,
        mean_expression = 0,
        median_expression = 0,
        max_expression = 0,
        n_expressing_cells = 0,
        total_cells = 0,
        pct_expressing = 0,
        expr_variance = 0,
        expr_cv = 0,
        expr_q10 = 0,
        expr_q25 = 0,
        expr_q50 = 0,
        expr_q75 = 0,
        expr_q90 = 0,
        expr_q95 = 0,
        high_expr_cells = 0,
        low_expr_cells = 0,
        pct_high_expr = 0,
        pct_low_expr = 0,
        stringsAsFactors = FALSE
      ))
    }
  }, error = function(e) {
    return(data.frame(
      sample = sample_name,
      donor_id = sample_name,
      mean_expression = 0,
      median_expression = 0,
      max_expression = 0,
      n_expressing_cells = 0,
      total_cells = 0,
      pct_expressing = 0,
      expr_variance = 0,
      expr_cv = 0,
      expr_q10 = 0,
      expr_q25 = 0,
      expr_q50 = 0,
      expr_q75 = 0,
      expr_q90 = 0,
      expr_q95 = 0,
      high_expr_cells = 0,
      low_expr_cells = 0,
      pct_high_expr = 0,
      pct_low_expr = 0,
      stringsAsFactors = FALSE
    ))
  })
}

# Get sample directories
sample_dirs <- list.dirs("data/raw/extracted", recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[!grepl("^\\.", sample_dirs)]

cat("Analyzing APOE expression in", length(sample_dirs), "donors...\n")

# Extract APOE data from all samples
all_results <- lapply(sample_dirs, extract_apoe_comprehensive)
all_results <- do.call(rbind, all_results)

# Combine with metadata
if (!is.null(all_results) && nrow(all_results) > 0) {
  all_results <- merge(all_results, sample_info, by.x = "sample", by.y = "sample_id", all.x = TRUE)
  
  # Save comprehensive results
  write.csv(all_results, "results/APOE_analysis/APOE_comprehensive_results.csv", row.names = FALSE)
  
  # Create comprehensive summary
  cat("\n=== APOE EXPRESSION SUMMARY ===\n")
  cat("Total donors analyzed:", length(unique(all_results$donor_id)), "\n")
  cat("Note: Each donor represents one sample (single-donor design)\n\n")
  
  # Overall summary
  overall_summary <- all_results %>%
    summarise(
      n_donors = n(),
      mean_expression = mean(mean_expression, na.rm = TRUE),
      median_expression = median(median_expression, na.rm = TRUE),
      mean_pct_expressing = mean(pct_expressing, na.rm = TRUE),
      mean_expr_cv = mean(expr_cv, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("Overall APOE Expression Summary:\n")
  print(overall_summary)
  
  # IDH status summary (similar to genotype stratification)
  idh_summary <- all_results %>%
    group_by(idh_status) %>%
    summarise(
      n_donors = n(),
      mean_expression = mean(mean_expression, na.rm = TRUE),
      median_expression = median(median_expression, na.rm = TRUE),
      mean_pct_expressing = mean(pct_expressing, na.rm = TRUE),
      mean_expr_cv = mean(expr_cv, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\nAPOE Expression by IDH Status:\n")
  print(idh_summary)
  
  # Tumor stage summary
  stage_summary <- all_results %>%
    group_by(tumor_stage) %>%
    summarise(
      n_donors = n(),
      mean_expression = mean(mean_expression, na.rm = TRUE),
      mean_pct_expressing = mean(pct_expressing, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\nAPOE Expression by Tumor Stage:\n")
  print(stage_summary)
  
  # Diagnosis summary
  diagnosis_summary <- all_results %>%
    group_by(diagnosis) %>%
    summarise(
      n_donors = n(),
      mean_expression = mean(mean_expression, na.rm = TRUE),
      mean_pct_expressing = mean(pct_expressing, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\nAPOE Expression by Diagnosis:\n")
  print(diagnosis_summary)
  
  # Create comprehensive visualizations
  if (nrow(all_results) > 0) {
    # APOE expression by IDH status
    p1 <- ggplot(all_results, aes(x = idh_status, y = mean_expression, fill = idh_status)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.6) +
      theme_minimal() +
      labs(title = "APOE Expression by IDH Status",
           subtitle = "Single-donor scRNA-seq samples",
           x = "IDH Status", y = "Mean Expression", fill = "IDH Status")
    
    ggsave("figures/APOE_analysis/APOE_by_IDH_status.pdf", p1, width = 10, height = 6)
    
    # APOE expression by tumor stage
    p2 <- ggplot(all_results, aes(x = tumor_stage, y = mean_expression, fill = tumor_stage)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.6) +
      theme_minimal() +
      labs(title = "APOE Expression by Tumor Stage",
           subtitle = "Single-donor scRNA-seq samples",
           x = "Tumor Stage", y = "Mean Expression", fill = "Stage")
    
    ggsave("figures/APOE_analysis/APOE_by_tumor_stage.pdf", p2, width = 10, height = 6)
    
    # APOE expression percentage by condition
    p3 <- ggplot(all_results, aes(x = condition_group, y = pct_expressing, fill = idh_status)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.6) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "APOE Expression Percentage by Condition",
           subtitle = "Percentage of cells expressing APOE",
           x = "Condition Group", y = "Percentage of Expressing Cells", fill = "IDH Status")
    
    ggsave("figures/APOE_analysis/APOE_percentage_by_condition.pdf", p3, width = 12, height = 6)
    
    # Expression distribution by diagnosis
    p4 <- ggplot(all_results, aes(x = diagnosis, y = mean_expression, fill = idh_status)) +
      geom_boxplot() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "APOE Expression by Diagnosis and IDH Status",
           subtitle = "Single-donor scRNA-seq samples",
           x = "Diagnosis", y = "Mean Expression", fill = "IDH Status")
    
    ggsave("figures/APOE_analysis/APOE_by_diagnosis.pdf", p4, width = 12, height = 6)
    
    # Expression coefficient of variation
    p5 <- ggplot(all_results, aes(x = idh_status, y = expr_cv, fill = idh_status)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.6) +
      theme_minimal() +
      labs(title = "APOE Expression Variability by IDH Status",
           subtitle = "Coefficient of Variation (CV)",
           x = "IDH Status", y = "Expression CV", fill = "IDH Status")
    
    ggsave("figures/APOE_analysis/APOE_variability_by_IDH.pdf", p5, width = 10, height = 6)
  }
  
  # Create summary report
  cat("\n=== ANALYSIS SUMMARY ===\n")
  cat("✓ APOE expression analyzed across", length(unique(all_results$donor_id)), "donors\n")
  cat("✓ Expression patterns compared by IDH status, tumor stage, and diagnosis\n")
  cat("✓ Comprehensive metrics calculated (mean, median, CV, percentiles)\n")
  cat("✓ Results saved to: results/APOE_analysis/APOE_comprehensive_results.csv\n")
  cat("✓ Figures saved to: figures/APOE_analysis/\n\n")
  
  cat("=== LIMITATIONS & RECOMMENDATIONS ===\n")
  cat("⚠️  No germline APOE genotype data (rs429358, rs7412) available\n")
  cat("⚠️  Cannot determine ε2/ε3/ε4 isoforms from RNA-seq alone\n")
  cat("💡 Recommendation: Obtain donor DNA genotypes for rs429358 and rs7412\n")
  cat("💡 Then re-analyze APOE expression by true genotype groups\n")
  
} else {
  cat("No results obtained from APOE analysis.\n")
}

cat("\n=== APOE BEST PRACTICE ANALYSIS COMPLETE ===\n")
