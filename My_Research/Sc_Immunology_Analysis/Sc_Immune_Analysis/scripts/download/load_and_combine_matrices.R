#!/usr/bin/env Rscript

# =============================================================================
# Load and Combine GSE222520 Expression Matrices
# Load 10x Genomics format data and combine into a single dataset
# =============================================================================

# Load required libraries
if (!requireNamespace("Matrix", quietly = TRUE))
    install.packages("Matrix")

if (!requireNamespace("readr", quietly = TRUE))
    install.packages("readr")

library(Matrix)
library(readr)
library(dplyr)

# Function to save script metadata
save_script_info <- function(script_name, description, output_files) {
    script_info <- list(
        script_name = script_name,
        description = description,
        timestamp = Sys.time(),
        output_files = output_files
    )
    
    info_file <- paste0("docs/data_description/", script_name, "_info.rds")
    saveRDS(script_info, info_file)
    cat("Saved script info to:", info_file, "\n")
}

# Function to read 10x Genomics data
read_10x_data <- function(data_path) {
    cat("Reading 10x data from:", data_path, "\n")
    
    # Read matrix
    matrix_file <- file.path(data_path, "matrix.mtx.gz")
    if (!file.exists(matrix_file)) {
        stop("Matrix file not found:", matrix_file)
    }
    
    # Read features (genes)
    features_file <- file.path(data_path, "features.tsv.gz")
    if (!file.exists(features_file)) {
        stop("Features file not found:", features_file)
    }
    
    # Read barcodes (cells)
    barcodes_file <- file.path(data_path, "barcodes.tsv.gz")
    if (!file.exists(barcodes_file)) {
        stop("Barcodes file not found:", barcodes_file)
    }
    
    # Read the data
    matrix_data <- Matrix::readMM(gzfile(matrix_file))
    features <- read_tsv(features_file, col_names = FALSE, show_col_types = FALSE)
    barcodes <- read_tsv(barcodes_file, col_names = FALSE, show_col_types = FALSE)
    
    # Set row and column names
    rownames(matrix_data) <- features$X1
    colnames(matrix_data) <- barcodes$X1
    
    # Return list with all components
    return(list(
        matrix = matrix_data,
        features = features,
        barcodes = barcodes
    ))
}

# Get list of extracted samples
extracted_samples <- list.dirs("data/raw/extracted", full.names = FALSE, recursive = FALSE)
extracted_samples <- extracted_samples[extracted_samples != ""]

cat("Found", length(extracted_samples), "extracted samples\n")

# Load each sample and collect data
all_data <- list()
sample_info <- data.frame()

for (i in seq_along(extracted_samples)) {
    sample <- extracted_samples[i]
    sample_dir <- file.path("data/raw/extracted", sample)
    
    cat("\nProcessing sample", i, "of", length(extracted_samples), ":", sample, "\n")
    
    # Check if the sample has the expected structure
    matrix_path <- file.path(sample_dir, "filtered_feature_bc_matrix")
    if (!dir.exists(matrix_path)) {
        cat("  Skipping - no filtered_feature_bc_matrix directory\n")
        next
    }
    
    tryCatch({
        # Read the 10x data
        data <- read_10x_data(matrix_path)
        
        # Add sample prefix to cell barcodes to avoid conflicts
        colnames(data$matrix) <- paste0(sample, "_", colnames(data$matrix))
        
        # Store the data
        all_data[[sample]] <- data
        
        # Record sample information
        sample_info <- rbind(sample_info, data.frame(
            sample = sample,
            sample_type = gsub(".*_", "", sample),
            n_cells = ncol(data$matrix),
            n_genes = nrow(data$matrix),
            total_counts = sum(data$matrix),
            mean_counts_per_cell = mean(colSums(data$matrix)),
            median_counts_per_cell = median(colSums(data$matrix))
        ))
        
        cat("  Loaded:", ncol(data$matrix), "cells,", nrow(data$matrix), "genes\n")
        
    }, error = function(e) {
        cat("  Error loading sample:", e$message, "\n")
    })
}

cat("\n=== LOADING SUMMARY ===\n")
cat("Successfully loaded", length(all_data), "samples\n")

if (nrow(sample_info) > 0) {
    print(sample_info)
    
    # Save sample information
    write.csv(sample_info, "data/processed/sample_info.csv", row.names = FALSE)
    saveRDS(sample_info, "data/processed/sample_info.rds")
}

# Combine all matrices if we have data
if (length(all_data) > 0) {
    cat("\n=== COMBINING MATRICES ===\n")
    
    # Get common genes across all samples
    all_genes <- lapply(all_data, function(x) rownames(x$matrix))
    common_genes <- Reduce(intersect, all_genes)
    cat("Common genes across all samples:", length(common_genes), "\n")
    
    # Combine matrices
    combined_matrix <- Matrix::Matrix(0, nrow = length(common_genes), ncol = 0, sparse = TRUE)
    rownames(combined_matrix) <- common_genes
    
    for (sample in names(all_data)) {
        cat("Adding sample:", sample, "\n")
        sample_matrix <- all_data[[sample]]$matrix[common_genes, , drop = FALSE]
        combined_matrix <- cbind(combined_matrix, sample_matrix)
    }
    
    cat("Combined matrix dimensions:", dim(combined_matrix), "\n")
    
    # Save combined matrix
    saveRDS(combined_matrix, "data/processed/GSE222520_combined_expression_matrix.rds")
    
    # Also save as a more accessible format (first 1000 genes for testing)
    if (nrow(combined_matrix) > 1000) {
        test_matrix <- combined_matrix[1:1000, ]
        write.csv(as.matrix(test_matrix), "data/processed/GSE222520_combined_expression_matrix_test.csv")
        cat("Saved test matrix (first 1000 genes) as CSV\n")
    }
    
    # Create cell metadata
    cell_metadata <- data.frame(
        cell_barcode = colnames(combined_matrix),
        sample = sapply(strsplit(colnames(combined_matrix), "_"), function(x) paste(x[1:(length(x)-1)], collapse = "_")),
        sample_type = sapply(strsplit(colnames(combined_matrix), "_"), function(x) x[length(x)-1])
    )
    
    # Add counts per cell
    cell_metadata$total_counts <- colSums(combined_matrix)
    cell_metadata$n_genes <- colSums(combined_matrix > 0)
    
    # Save cell metadata
    write.csv(cell_metadata, "data/processed/GSE222520_cell_metadata.csv", row.names = FALSE)
    saveRDS(cell_metadata, "data/processed/GSE222520_cell_metadata.rds")
    
    # Create gene metadata
    gene_metadata <- data.frame(
        gene_symbol = rownames(combined_matrix),
        total_counts = rowSums(combined_matrix),
        n_cells_expressed = rowSums(combined_matrix > 0),
        mean_expression = rowMeans(combined_matrix)
    )
    
    # Save gene metadata
    write.csv(gene_metadata, "data/processed/GSE222520_gene_metadata.csv", row.names = FALSE)
    saveRDS(gene_metadata, "data/processed/GSE222520_gene_metadata.rds")
    
    cat("\n=== COMBINED DATASET SUMMARY ===\n")
    cat("Total cells:", ncol(combined_matrix), "\n")
    cat("Total genes:", nrow(combined_matrix), "\n")
    cat("Total UMIs:", sum(combined_matrix), "\n")
    cat("Mean UMIs per cell:", mean(colSums(combined_matrix)), "\n")
    cat("Median UMIs per cell:", median(colSums(combined_matrix)), "\n")
    
    # Sample type distribution
    cat("\nSample type distribution:\n")
    print(table(cell_metadata$sample_type))
    
} else {
    cat("No data to combine\n")
}

# Save script information
save_script_info(
    script_name = "load_and_combine_matrices.R",
    description = "Load and combine 10x Genomics expression matrices",
    output_files = list.files("data/processed", recursive = TRUE, full.names = TRUE)
)

cat("\n=== LOADING AND COMBINING COMPLETE ===\n")
cat("Check data/processed/ for combined datasets\n")
