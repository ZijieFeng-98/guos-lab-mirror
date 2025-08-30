# Simple APOE Isoform Analysis in Single-Cell Data
# GSE222520 - Brain Tumor Leukocyte Single-Cell Analysis
# Author: Data Analysis Pipeline
# Date: 2024

# Load required libraries (basic R packages)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(Matrix)

# Set random seed for reproducibility
set.seed(42)

# Create output directories
dir.create("results/APOE_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/APOE_analysis", recursive = TRUE, showWarnings = FALSE)

# Function to read 10x data without Seurat
read_10x_data <- function(data_path, sample_name) {
  cat("Reading data for sample:", sample_name, "\n")
  
  # Check if filtered_feature_bc_matrix exists
  matrix_path <- file.path(data_path, "filtered_feature_bc_matrix")
  if (!dir.exists(matrix_path)) {
    # Look for the actual matrix directory
    subdirs <- list.dirs(data_path, recursive = FALSE)
    matrix_path <- file.path(subdirs[1], "filtered_feature_bc_matrix")
  }
  
  if (!dir.exists(matrix_path)) {
    stop(paste("Matrix directory not found for sample:", sample_name))
  }
  
  # Read the matrix files
  barcodes_file <- file.path(matrix_path, "barcodes.tsv.gz")
  features_file <- file.path(matrix_path, "features.tsv.gz")
  matrix_file <- file.path(matrix_path, "matrix.mtx.gz")
  
  # Read barcodes (cell IDs)
  barcodes <- readLines(gzfile(barcodes_file))
  
  # Read features (gene names)
  features <- readLines(gzfile(features_file))
  gene_names <- sapply(strsplit(features, "\t"), function(x) x[2])
  
  # Read matrix
  matrix_data <- readLines(gzfile(matrix_file))
  
  # Parse matrix header
  header <- strsplit(matrix_data[1], " ")[[1]]
  n_genes <- as.numeric(header[1])
  n_cells <- as.numeric(header[2])
  n_entries <- as.numeric(header[3])
  
  cat("  Matrix dimensions:", n_genes, "genes x", n_cells, "cells\n")
  
  # Parse matrix entries
  entries <- matrix_data[2:(n_entries + 1)]
  parsed_entries <- do.call(rbind, lapply(entries, function(x) {
    parts <- strsplit(x, " ")[[1]]
    c(as.numeric(parts[1]), as.numeric(parts[2]), as.numeric(parts[3]))
  }))
  
  # Create sparse matrix
  sparse_matrix <- sparseMatrix(
    i = parsed_entries[, 1],
    j = parsed_entries[, 2],
    x = parsed_entries[, 3],
    dims = c(n_genes, n_cells),
    dimnames = list(gene_names, barcodes)
  )
  
  return(list(
    matrix = sparse_matrix,
    genes = gene_names,
    cells = barcodes,
    sample = sample_name
  ))
}

# Function to extract sample information from metadata
extract_sample_info <- function(metadata_file) {
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  
  # Extract key information
  sample_info <- data.frame(
    sample_id = metadata$title,
    geo_accession = metadata$geo_accession,
    tissue = metadata$`tissue:ch1`,
    diagnosis = metadata$`diagnosis:ch1`,
    genotype = metadata$`genotype:ch1`,
    primary_recurrent = metadata$`primary/recurrent:ch1`,
    treatment = metadata$`treatment:ch1`,
    cell_type = metadata$`cell type:ch1`,
    stringsAsFactors = FALSE
  )
  
  # Clean up sample names
  sample_info$sample_id <- gsub(", replicate1, scRNAseq", "", sample_info$sample_id)
  
  return(sample_info)
}

# Load sample metadata
cat("Loading sample metadata...\n")
sample_info <- extract_sample_info("data/metadata/GSE222520_phenotype_data.csv")

# Define sample groups
sample_info$group <- case_when(
  grepl("^NGB", sample_info$sample_id) ~ "Normal",
  grepl("^IMP", sample_info$sample_id) ~ "Primary_IDH_Mutant",
  grepl("^IWP", sample_info$sample_id) ~ "Primary_IDH_WildType",
  grepl("^IMR", sample_info$sample_id) ~ "Recurrent_IDH_Mutant",
  grepl("^IWR", sample_info$sample_id) ~ "Recurrent_IDH_WildType"
)

# List of sample directories
sample_dirs <- list.dirs("data/raw/extracted", recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[!grepl("^\\.", sample_dirs)]  # Remove hidden files

cat("Found", length(sample_dirs), "sample directories\n")

# Load all samples
sample_data <- list()

for (sample_dir in sample_dirs) {
  sample_name <- gsub("GSM[0-9]+_", "", sample_dir)
  sample_path <- file.path("data/raw/extracted", sample_dir)
  
  tryCatch({
    data <- read_10x_data(sample_path, sample_name)
    sample_data[[sample_name]] <- data
    cat("Successfully loaded:", sample_name, "with", ncol(data$matrix), "cells\n")
  }, error = function(e) {
    cat("Error loading", sample_name, ":", e$message, "\n")
  })
}

# APOE Analysis
cat("\n=== APOE ISOFORM ANALYSIS ===\n")

# Define APOE-related genes to look for
apoe_genes <- c("APOE", "APOE2", "APOE3", "APOE4", "APOC1", "APOC2", "APOC3", "APOC4")

# Function to analyze APOE expression in a sample
analyze_apoe_in_sample <- function(data, sample_name) {
  cat("Analyzing APOE in sample:", sample_name, "\n")
  
  # Get gene names
  genes <- data$genes
  
  # Find APOE-related genes
  apoe_found <- genes[grepl("APO[EC]", genes, ignore.case = TRUE)]
  
  cat("  APOE-related genes found:", paste(apoe_found, collapse = ", "), "\n")
  
  if (length(apoe_found) > 0) {
    # Calculate expression statistics
    matrix <- data$matrix
    
    apoe_stats <- data.frame(
      gene = apoe_found,
      sample = sample_name,
      total_cells = ncol(matrix),
      expressing_cells = sapply(apoe_found, function(gene) {
        gene_idx <- which(genes == gene)
        if (length(gene_idx) > 0) {
          sum(matrix[gene_idx, ] > 0)
        } else {
          0
        }
      }),
      mean_expression = sapply(apoe_found, function(gene) {
        gene_idx <- which(genes == gene)
        if (length(gene_idx) > 0) {
          mean(as.numeric(matrix[gene_idx, ]))
        } else {
          0
        }
      }),
      median_expression = sapply(apoe_found, function(gene) {
        gene_idx <- which(genes == gene)
        if (length(gene_idx) > 0) {
          median(as.numeric(matrix[gene_idx, ]))
        } else {
          0
        }
      }),
      max_expression = sapply(apoe_found, function(gene) {
        gene_idx <- which(genes == gene)
        if (length(gene_idx) > 0) {
          max(as.numeric(matrix[gene_idx, ]))
        } else {
          0
        }
      }),
      stringsAsFactors = FALSE
    )
    
    apoe_stats$percent_expressing <- (apoe_stats$expressing_cells / apoe_stats$total_cells) * 100
    
    return(apoe_stats)
  } else {
    return(NULL)
  }
}

# Analyze APOE in all samples
apoe_results <- list()
for (sample_name in names(sample_data)) {
  result <- analyze_apoe_in_sample(sample_data[[sample_name]], sample_name)
  if (!is.null(result)) {
    apoe_results[[sample_name]] <- result
  }
}

# Combine all results
if (length(apoe_results) > 0) {
  all_apoe_stats <- do.call(rbind, apoe_results)
  
  # Add sample metadata
  sample_meta <- sample_info[match(all_apoe_stats$sample, sample_info$sample_id), ]
  all_apoe_stats$group <- sample_meta$group
  all_apoe_stats$diagnosis <- sample_meta$diagnosis
  all_apoe_stats$genotype <- sample_meta$genotype
  all_apoe_stats$primary_recurrent <- sample_meta$primary_recurrent
  
  # Save APOE statistics
  write.csv(all_apoe_stats, "results/APOE_analysis/APOE_expression_stats.csv", row.names = FALSE)
  
  cat("\nAPOE expression statistics saved to: results/APOE_analysis/APOE_expression_stats.csv\n")
  
  # Create visualizations
  cat("Creating APOE expression visualizations...\n")
  
  # 1. APOE expression by sample
  if ("APOE" %in% all_apoe_stats$gene) {
    apoe_data <- all_apoe_stats[all_apoe_stats$gene == "APOE", ]
    
    p1 <- ggplot(apoe_data, aes(x = sample, y = mean_expression, fill = group)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "APOE Mean Expression by Sample",
           x = "Sample", y = "Mean Expression", fill = "Group")
    
    ggsave("figures/APOE_analysis/APOE_mean_expression_by_sample.pdf", p1, width = 12, height = 6)
    
    # 2. APOE expression by group
    p2 <- ggplot(apoe_data, aes(x = group, y = mean_expression, fill = group)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.6) +
      theme_minimal() +
      labs(title = "APOE Expression by Group",
           x = "Group", y = "Mean Expression", fill = "Group")
    
    ggsave("figures/APOE_analysis/APOE_expression_by_group.pdf", p2, width = 10, height = 6)
    
    # 3. APOE expression by diagnosis
    p3 <- ggplot(apoe_data, aes(x = diagnosis, y = mean_expression, fill = diagnosis)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.6) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "APOE Expression by Diagnosis",
           x = "Diagnosis", y = "Mean Expression", fill = "Diagnosis")
    
    ggsave("figures/APOE_analysis/APOE_expression_by_diagnosis.pdf", p3, width = 10, height = 6)
    
    # 4. APOE expression by genotype
    p4 <- ggplot(apoe_data, aes(x = genotype, y = mean_expression, fill = genotype)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.6) +
      theme_minimal() +
      labs(title = "APOE Expression by Genotype",
           x = "Genotype", y = "Mean Expression", fill = "Genotype")
    
    ggsave("figures/APOE_analysis/APOE_expression_by_genotype.pdf", p4, width = 8, height = 6)
    
    # 5. Percentage of cells expressing APOE
    p5 <- ggplot(apoe_data, aes(x = sample, y = percent_expressing, fill = group)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Percentage of Cells Expressing APOE",
           x = "Sample", y = "Percentage of Cells (%)", fill = "Group")
    
    ggsave("figures/APOE_analysis/APOE_percent_expressing_by_sample.pdf", p5, width = 12, height = 6)
  }
  
} else {
  cat("No APOE-related genes found in any sample.\n")
}

# Detailed cell-level analysis for samples with APOE expression
cat("\n=== DETAILED CELL-LEVEL APOE ANALYSIS ===\n")

# Function to analyze APOE expression at cell level
analyze_apoe_cell_level <- function(data, sample_name) {
  cat("Performing cell-level APOE analysis for:", sample_name, "\n")
  
  # Check if APOE is present
  if (!"APOE" %in% data$genes) {
    cat("  APOE not found in this sample\n")
    return(NULL)
  }
  
  # Get APOE expression
  gene_idx <- which(data$genes == "APOE")
  apoe_expr <- as.numeric(data$matrix[gene_idx, ])
  
  # Calculate cell-level metrics
  n_genes_per_cell <- colSums(data$matrix > 0)
  n_umis_per_cell <- colSums(data$matrix)
  
  # Create cell-level data
  cell_data <- data.frame(
    sample = sample_name,
    cell_id = data$cells,
    apoe_expression = apoe_expr,
    apoe_detected = apoe_expr > 0,
    n_genes = n_genes_per_cell,
    n_umis = n_umis_per_cell,
    stringsAsFactors = FALSE
  )
  
  return(cell_data)
}

# Perform cell-level analysis
cell_level_results <- list()
for (sample_name in names(sample_data)) {
  result <- analyze_apoe_cell_level(sample_data[[sample_name]], sample_name)
  if (!is.null(result)) {
    cell_level_results[[sample_name]] <- result
  }
}

# Combine cell-level results
if (length(cell_level_results) > 0) {
  all_cell_data <- do.call(rbind, cell_level_results)
  
  # Add sample metadata
  sample_meta <- sample_info[match(all_cell_data$sample, sample_info$sample_id), ]
  all_cell_data$group <- sample_meta$group
  all_cell_data$diagnosis <- sample_meta$diagnosis
  all_cell_data$genotype <- sample_meta$genotype
  
  # Save cell-level data
  write.csv(all_cell_data, "results/APOE_analysis/APOE_cell_level_data.csv", row.names = FALSE)
  
  # Create cell-level visualizations
  cat("Creating cell-level APOE visualizations...\n")
  
  # 1. APOE expression distribution by group
  p1 <- ggplot(all_cell_data[all_cell_data$apoe_expression > 0, ], 
               aes(x = group, y = log2(apoe_expression + 1), fill = group)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "APOE Expression Distribution (Expressing Cells Only)",
         x = "Group", y = "Log2(APOE Expression + 1)", fill = "Group")
  
  ggsave("figures/APOE_analysis/APOE_expression_distribution_by_group.pdf", p1, width = 10, height = 6)
  
  # 2. APOE expression vs number of genes
  p2 <- ggplot(all_cell_data[all_cell_data$apoe_expression > 0, ], 
               aes(x = n_genes, y = apoe_expression, color = group)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_minimal() +
    labs(title = "APOE Expression vs Number of Genes",
         x = "Number of Genes", y = "APOE Expression", color = "Group")
  
  ggsave("figures/APOE_analysis/APOE_vs_genes_scatter.pdf", p2, width = 10, height = 6)
  
  # 3. Percentage of APOE-expressing cells by group
  group_summary <- all_cell_data %>%
    group_by(group) %>%
    summarise(
      total_cells = n(),
      apoe_positive = sum(apoe_detected),
      percent_positive = (apoe_positive / total_cells) * 100,
      mean_expression = mean(apoe_expression),
      median_expression = median(apoe_expression)
    )
  
  write.csv(group_summary, "results/APOE_analysis/APOE_group_summary.csv", row.names = FALSE)
  
  p3 <- ggplot(group_summary, aes(x = group, y = percent_positive, fill = group)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Percentage of APOE-Expressing Cells by Group",
         x = "Group", y = "Percentage of Cells (%)", fill = "Group")
  
  ggsave("figures/APOE_analysis/APOE_percent_positive_by_group.pdf", p3, width = 10, height = 6)
  
  cat("Cell-level analysis complete. Results saved to results/APOE_analysis/\n")
}

# Statistical analysis
cat("\n=== STATISTICAL ANALYSIS ===\n")

if (length(apoe_results) > 0 && "APOE" %in% all_apoe_stats$gene) {
  apoe_data <- all_apoe_stats[all_apoe_stats$gene == "APOE", ]
  
  # ANOVA for group differences
  if (length(unique(apoe_data$group)) > 1) {
    aov_result <- aov(mean_expression ~ group, data = apoe_data)
    cat("ANOVA results for APOE expression by group:\n")
    print(summary(aov_result))
    
    # Save ANOVA results
    sink("results/APOE_analysis/APOE_anova_results.txt")
    cat("ANOVA results for APOE expression by group:\n")
    print(summary(aov_result))
    sink()
  }
  
  # T-test for specific comparisons
  if ("Normal" %in% apoe_data$group) {
    normal_data <- apoe_data$mean_expression[apoe_data$group == "Normal"]
    
    for (group in unique(apoe_data$group[apoe_data$group != "Normal"])) {
      group_data <- apoe_data$mean_expression[apoe_data$group == group]
      
      if (length(normal_data) > 0 && length(group_data) > 0) {
        t_result <- t.test(normal_data, group_data)
        cat("\nT-test: Normal vs", group, "\n")
        print(t_result)
      }
    }
  }
}

# Summary
cat("\n=== APOE ANALYSIS SUMMARY ===\n")
cat("Total samples analyzed:", length(sample_data), "\n")

if (length(apoe_results) > 0) {
  cat("Samples with APOE expression:", length(apoe_results), "\n")
  cat("APOE-related genes found:", length(unique(all_apoe_stats$gene)), "\n")
  cat("Results saved to: results/APOE_analysis/\n")
  cat("Figures saved to: figures/APOE_analysis/\n")
} else {
  cat("No APOE expression detected in any sample.\n")
}

cat("\n=== APOE ANALYSIS COMPLETE ===\n")
