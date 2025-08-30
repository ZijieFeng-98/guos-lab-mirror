# APOE Isoform Analysis in Single-Cell Data
# GSE222520 - Brain Tumor Leukocyte Single-Cell Analysis
# Author: Data Analysis Pipeline
# Date: 2024

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(tidyr)
library(stringr)
library(Matrix)

# Set random seed for reproducibility
set.seed(42)

# Create output directories
dir.create("results/APOE_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/APOE_analysis", recursive = TRUE, showWarnings = FALSE)

# Function to load 10x data
load_10x_data <- function(data_path, sample_name) {
  cat("Loading data for sample:", sample_name, "\n")
  
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
  
  # Load the data
  data <- Read10X(matrix_path)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = sample_name,
    min.cells = 1,  # Include all genes for APOE analysis
    min.features = 1  # Include all cells for APOE analysis
  )
  
  # Add sample metadata
  seurat_obj$sample_id <- sample_name
  
  return(seurat_obj)
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
seurat_objects <- list()
sample_paths <- list()

for (sample_dir in sample_dirs) {
  sample_name <- gsub("GSM[0-9]+_", "", sample_dir)
  sample_path <- file.path("data/raw/extracted", sample_dir)
  
  tryCatch({
    seurat_obj <- load_10x_data(sample_path, sample_name)
    seurat_objects[[sample_name]] <- seurat_obj
    sample_paths[[sample_name]] <- sample_path
    cat("Successfully loaded:", sample_name, "with", ncol(seurat_obj), "cells\n")
  }, error = function(e) {
    cat("Error loading", sample_name, ":", e$message, "\n")
  })
}

# APOE Analysis
cat("\n=== APOE ISOFORM ANALYSIS ===\n")

# Define APOE-related genes to look for
apoe_genes <- c("APOE", "APOE2", "APOE3", "APOE4", "APOC1", "APOC2", "APOC3", "APOC4")

# Function to analyze APOE expression in a sample
analyze_apoe_in_sample <- function(seurat_obj, sample_name) {
  cat("Analyzing APOE in sample:", sample_name, "\n")
  
  # Get gene names
  genes <- rownames(seurat_obj)
  
  # Find APOE-related genes
  apoe_found <- genes[grepl("APO[EC]", genes, ignore.case = TRUE)]
  
  cat("  APOE-related genes found:", paste(apoe_found, collapse = ", "), "\n")
  
  if (length(apoe_found) > 0) {
    # Calculate expression statistics
    expr_data <- GetAssayData(seurat_obj, slot = "counts")
    
    apoe_stats <- data.frame(
      gene = apoe_found,
      sample = sample_name,
      total_cells = ncol(seurat_obj),
      expressing_cells = sapply(apoe_found, function(gene) {
        sum(expr_data[gene, ] > 0)
      }),
      mean_expression = sapply(apoe_found, function(gene) {
        mean(expr_data[gene, ])
      }),
      median_expression = sapply(apoe_found, function(gene) {
        median(expr_data[gene, ])
      }),
      max_expression = sapply(apoe_found, function(gene) {
        max(expr_data[gene, ])
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
for (sample_name in names(seurat_objects)) {
  result <- analyze_apoe_in_sample(seurat_objects[[sample_name]], sample_name)
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
  
  # 6. All APOE-related genes heatmap
  if (nrow(all_apoe_stats) > 1) {
    # Create a heatmap of all APOE-related genes
    heatmap_data <- all_apoe_stats %>%
      select(gene, sample, mean_expression, group) %>%
      spread(sample, mean_expression, fill = 0)
    
    # Convert to matrix for heatmap
    heatmap_matrix <- as.matrix(heatmap_data[, -c(1, 2)])
    rownames(heatmap_matrix) <- heatmap_data$gene
    
    # Create heatmap
    pdf("figures/APOE_analysis/APOE_genes_heatmap.pdf", width = 12, height = 8)
    pheatmap::pheatmap(heatmap_matrix, 
                       scale = "row",
                       annotation_col = data.frame(
                         Group = heatmap_data$group,
                         row.names = colnames(heatmap_matrix)
                       ),
                       main = "APOE-related Genes Expression Heatmap",
                       fontsize_row = 10,
                       fontsize_col = 8)
    dev.off()
  }
  
} else {
  cat("No APOE-related genes found in any sample.\n")
}

# Detailed cell-level analysis for samples with APOE expression
cat("\n=== DETAILED CELL-LEVEL APOE ANALYSIS ===\n")

# Function to analyze APOE expression at cell level
analyze_apoe_cell_level <- function(seurat_obj, sample_name) {
  cat("Performing cell-level APOE analysis for:", sample_name, "\n")
  
  # Check if APOE is present
  if (!"APOE" %in% rownames(seurat_obj)) {
    cat("  APOE not found in this sample\n")
    return(NULL)
  }
  
  # Get APOE expression
  apoe_expr <- GetAssayData(seurat_obj, slot = "counts")["APOE", ]
  
  # Create cell-level data
  cell_data <- data.frame(
    sample = sample_name,
    cell_id = colnames(seurat_obj),
    apoe_expression = as.numeric(apoe_expr),
    apoe_detected = as.numeric(apoe_expr) > 0,
    n_genes = seurat_obj$nFeature_RNA,
    n_umis = seurat_obj$nCount_RNA,
    stringsAsFactors = FALSE
  )
  
  return(cell_data)
}

# Perform cell-level analysis
cell_level_results <- list()
for (sample_name in names(seurat_objects)) {
  result <- analyze_apoe_cell_level(seurat_objects[[sample_name]], sample_name)
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
cat("Total samples analyzed:", length(seurat_objects), "\n")

if (length(apoe_results) > 0) {
  cat("Samples with APOE expression:", length(apoe_results), "\n")
  cat("APOE-related genes found:", length(unique(all_apoe_stats$gene)), "\n")
  cat("Results saved to: results/APOE_analysis/\n")
  cat("Figures saved to: figures/APOE_analysis/\n")
} else {
  cat("No APOE expression detected in any sample.\n")
}

cat("\n=== APOE ANALYSIS COMPLETE ===\n")
