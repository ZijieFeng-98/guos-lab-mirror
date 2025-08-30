#!/usr/bin/env Rscript

# Download and Organize GSE222520 Dataset
# Single-cell RNA-seq analysis of immune cells

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE))
    BiocManager::install("GEOquery")

if (!requireNamespace("Biobase", quietly = TRUE))
    BiocManager::install("Biobase")

library(GEOquery)
library(Biobase)

# Create project directories
dirs <- c("data", "data/raw", "data/processed", "scripts", "results", "figures", "docs")
for (dir in dirs) {
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
        cat("Created directory:", dir, "\n")
    }
}

# Download GSE222520 dataset
cat("Downloading GSE222520 dataset from GEO...\n")
gse <- getGEO("GSE222520", GSEMatrix = TRUE, AnnotGPL = TRUE)

# Check what we downloaded
cat("Dataset structure:\n")
print(length(gse))
print(names(gse))

# Save the raw GEO object
saveRDS(gse, "data/raw/GSE222520_raw.rds")
cat("Saved raw GEO object to data/raw/GSE222520_raw.rds\n")

# Extract expression data and metadata
if (length(gse) > 0) {
    # Get the first dataset (usually there's only one)
    eset <- gse[[1]]
    
    # Extract expression matrix
    expr_matrix <- exprs(eset)
    cat("Expression matrix dimensions:", dim(expr_matrix), "\n")
    
    # Extract phenotype data (metadata)
    pheno_data <- pData(eset)
    cat("Phenotype data dimensions:", dim(pheno_data), "\n")
    
    # Extract feature data (gene annotations)
    feature_data <- fData(eset)
    cat("Feature data dimensions:", dim(feature_data), "\n")
    
    # Save processed data
    saveRDS(expr_matrix, "data/processed/GSE222520_expression_matrix.rds")
    saveRDS(pheno_data, "data/processed/GSE222520_phenotype_data.rds")
    saveRDS(feature_data, "data/processed/GSE222520_feature_data.rds")
    
    # Also save as CSV for easier inspection
    write.csv(expr_matrix, "data/processed/GSE222520_expression_matrix.csv")
    write.csv(pheno_data, "data/processed/GSE222520_phenotype_data.csv")
    write.csv(feature_data, "data/processed/GSE222520_feature_data.csv")
    
    cat("Saved processed data to data/processed/\n")
    
    # Print summary information
    cat("\n=== DATASET SUMMARY ===\n")
    cat("GEO Accession: GSE222520\n")
    cat("Expression matrix: ", nrow(expr_matrix), "genes x ", ncol(expr_matrix), "cells\n")
    cat("Phenotype variables: ", ncol(pheno_data), "\n")
    cat("Feature annotations: ", ncol(feature_data), "\n")
    
    # Show available phenotype variables
    cat("\nAvailable phenotype variables:\n")
    print(colnames(pheno_data))
    
} else {
    cat("No data found in GSE222520\n")
}

cat("\nDownload and organization complete!\n")
