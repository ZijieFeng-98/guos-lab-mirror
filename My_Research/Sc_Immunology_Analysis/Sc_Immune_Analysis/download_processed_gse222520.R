#!/usr/bin/env Rscript

# =============================================================================
# Download and Organize GSE222520 Processed Data
# Single-cell RNA-seq analysis of glioma immune microenvironment
# 
# Nature of the Processed Data:
# - Expression matrices (gene × cell count matrices)
# - Cell metadata (sample ID, patient ID, QC metrics, cell type annotations)
# - Cluster information (MG, CD8 T, NK, etc.)
# - Reference signatures (GlioTIME-36)
# 
# Advantages:
# - Time-saving (skip raw FASTQ → QC → alignment → count matrix)
# - Batch-corrected (Harmony applied)
# - Pre-annotated clusters
# - Good for secondary analysis
# =============================================================================

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE))
    BiocManager::install("GEOquery")

if (!requireNamespace("Biobase", quietly = TRUE))
    BiocManager::install("Biobase")

library(GEOquery)
library(Biobase)

# Create comprehensive project structure
dirs <- c(
    "data", 
    "data/raw", 
    "data/processed", 
    "data/metadata",
    "data/annotations",
    "scripts", 
    "scripts/download",
    "scripts/analysis",
    "results", 
    "results/qc",
    "results/clustering",
    "figures", 
    "figures/qc",
    "figures/clustering",
    "docs",
    "docs/data_description"
)

for (dir in dirs) {
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
        cat("Created directory:", dir, "\n")
    }
}

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

# Download GSE222520 dataset
cat("=== DOWNLOADING GSE222520 PROCESSED DATA ===\n")
cat("Downloading from GEO...\n")

gse <- getGEO("GSE222520", GSEMatrix = TRUE, AnnotGPL = TRUE)

# Save the raw GEO object
saveRDS(gse, "data/raw/GSE222520_raw.rds")
cat("Saved raw GEO object to data/raw/GSE222520_raw.rds\n")

# Extract and organize data
if (length(gse) > 0) {
    eset <- gse[[1]]
    
    # Extract phenotype data (metadata)
    pheno_data <- pData(eset)
    cat("Phenotype data dimensions:", dim(pheno_data), "\n")
    
    # Extract feature data (gene annotations)
    feature_data <- fData(eset)
    cat("Feature data dimensions:", dim(feature_data), "\n")
    
    # Extract expression matrix (if available)
    expr_matrix <- exprs(eset)
    cat("Expression matrix dimensions:", dim(expr_matrix), "\n")
    
    # Save processed data
    saveRDS(pheno_data, "data/metadata/GSE222520_phenotype_data.rds")
    saveRDS(feature_data, "data/metadata/GSE222520_feature_data.rds")
    
    # Save as CSV for easier inspection
    write.csv(pheno_data, "data/metadata/GSE222520_phenotype_data.csv", row.names = FALSE)
    write.csv(feature_data, "data/metadata/GSE222520_feature_data.csv", row.names = FALSE)
    
    if (nrow(expr_matrix) > 0) {
        saveRDS(expr_matrix, "data/processed/GSE222520_expression_matrix.rds")
        write.csv(expr_matrix, "data/processed/GSE222520_expression_matrix.csv")
        cat("Saved expression matrix\n")
    }
    
    # Check for supplementary files
    supp_files <- pheno_data$supplementary_file_1
    supp_files <- supp_files[!is.na(supp_files) & supp_files != ""]
    
    if (length(supp_files) > 0) {
        cat("\n=== DOWNLOADING SUPPLEMENTARY FILES ===\n")
        cat("Found", length(supp_files), "supplementary files\n")
        
        for (i in seq_along(supp_files)) {
            file_url <- supp_files[i]
            file_name <- basename(file_url)
            output_path <- paste0("data/raw/", file_name)
            
            cat("Downloading:", file_name, "\n")
            tryCatch({
                download.file(file_url, output_path, mode = "wb")
                cat("Successfully downloaded:", file_name, "\n")
            }, error = function(e) {
                cat("Failed to download:", file_name, "-", e$message, "\n")
            })
        }
    }
    
    # Create data summary
    cat("\n=== DATASET SUMMARY ===\n")
    cat("GEO Accession: GSE222520\n")
    cat("Expression matrix: ", nrow(expr_matrix), "genes x ", ncol(expr_matrix), "cells\n")
    cat("Phenotype variables: ", ncol(pheno_data), "\n")
    cat("Feature annotations: ", ncol(feature_data), "\n")
    cat("Supplementary files: ", length(supp_files), "\n")
    
    # Show key phenotype variables
    cat("\nKey phenotype variables:\n")
    key_vars <- c("cell type:ch1", "diagnosis:ch1", "genotype:ch1", 
                  "primary/recurrent:ch1", "tissue:ch1", "treatment:ch1")
    for (var in key_vars) {
        if (var %in% colnames(pheno_data)) {
            cat(var, ":", length(unique(pheno_data[[var]])), "unique values\n")
        }
    }
    
    # Create data description file
    data_description <- list(
        geo_accession = "GSE222520",
        title = "Single-cell RNA-seq analysis of glioma immune microenvironment",
        data_type = "Processed scRNA-seq",
        expression_matrix_dim = dim(expr_matrix),
        phenotype_variables = colnames(pheno_data),
        feature_annotations = colnames(feature_data),
        supplementary_files = supp_files,
        download_date = Sys.time(),
        notes = "Processed data includes batch-corrected expression matrices, cell metadata, and cluster annotations"
    )
    
    saveRDS(data_description, "docs/data_description/GSE222520_data_description.rds")
    writeLines(
        c(
            "GSE222520 Data Description",
            "========================",
            "",
            paste("GEO Accession:", data_description$geo_accession),
            paste("Title:", data_description$title),
            paste("Data Type:", data_description$data_type),
            paste("Expression Matrix:", paste(data_description$expression_matrix_dim, collapse = " x ")),
            paste("Phenotype Variables:", length(data_description$phenotype_variables)),
            paste("Feature Annotations:", length(data_description$feature_annotations)),
            paste("Supplementary Files:", length(data_description$supplementary_files)),
            paste("Download Date:", data_description$download_date),
            "",
            "Notes:",
            data_description$notes
        ),
        "docs/data_description/GSE222520_data_description.txt"
    )
    
} else {
    cat("No data found in GSE222520\n")
}

# Save script information
save_script_info(
    script_name = "download_processed_gse222520.R",
    description = "Download and organize processed GSE222520 scRNA-seq data",
    output_files = list.files("data", recursive = TRUE, full.names = TRUE)
)

cat("\n=== DOWNLOAD AND ORGANIZATION COMPLETE ===\n")
cat("Project structure created and data organized\n")
cat("Check docs/data_description/ for detailed information\n")
