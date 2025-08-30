#!/usr/bin/env Rscript

# =============================================================================
# Extract and Organize GSE222520 Processed Data
# Extract tar.gz files and organize expression matrices and metadata
# =============================================================================

# Load required libraries
library(dplyr)
library(readr)

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

# Create extraction directories
extract_dirs <- c(
    "data/raw/extracted",
    "data/processed/expression_matrices",
    "data/processed/cell_metadata",
    "data/processed/gene_metadata"
)

for (dir in extract_dirs) {
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
        cat("Created directory:", dir, "\n")
    }
}

# Get list of downloaded tar.gz files
tar_files <- list.files("data/raw", pattern = "\\.tar\\.gz$", full.names = TRUE)
cat("Found", length(tar_files), "tar.gz files to extract\n")

# Extract each tar.gz file
extracted_samples <- c()

for (tar_file in tar_files) {
    sample_name <- gsub("\\.tar\\.gz$", "", basename(tar_file))
    extract_dir <- paste0("data/raw/extracted/", sample_name)
    
    cat("Extracting:", basename(tar_file), "\n")
    
    # Create sample-specific directory
    if (!dir.exists(extract_dir)) {
        dir.create(extract_dir, recursive = TRUE)
    }
    
    # Extract tar.gz file
    tryCatch({
        untar(tar_file, exdir = extract_dir)
        extracted_samples <- c(extracted_samples, sample_name)
        cat("Successfully extracted:", sample_name, "\n")
    }, error = function(e) {
        cat("Failed to extract:", sample_name, "-", e$message, "\n")
    })
}

cat("\nExtracted", length(extracted_samples), "samples\n")

# Look for common file patterns in extracted data
cat("\n=== EXPLORING EXTRACTED DATA STRUCTURE ===\n")

# Check what files are in each extracted directory
for (sample in extracted_samples) {
    sample_dir <- paste0("data/raw/extracted/", sample)
    files <- list.files(sample_dir, recursive = TRUE, full.names = TRUE)
    
    cat("\nSample:", sample, "-", length(files), "files\n")
    
    # Look for common scRNA-seq file patterns
    mtx_files <- files[grepl("\\.mtx$", files)]
    csv_files <- files[grepl("\\.csv$", files)]
    tsv_files <- files[grepl("\\.tsv$", files)]
    h5_files <- files[grepl("\\.h5$", files)]
    rds_files <- files[grepl("\\.rds$", files)]
    
    if (length(mtx_files) > 0) cat("  MTX files:", length(mtx_files), "\n")
    if (length(csv_files) > 0) cat("  CSV files:", length(csv_files), "\n")
    if (length(tsv_files) > 0) cat("  TSV files:", length(tsv_files), "\n")
    if (length(h5_files) > 0) cat("  H5 files:", length(h5_files), "\n")
    if (length(rds_files) > 0) cat("  RDS files:", length(rds_files), "\n")
    
    # Show first few files
    if (length(files) > 0) {
        cat("  First few files:\n")
        for (i in 1:min(5, length(files))) {
            cat("    ", basename(files[i]), "\n")
        }
    }
}

# Create a summary of all extracted data
extraction_summary <- list(
    total_samples = length(extracted_samples),
    extracted_samples = extracted_samples,
    extraction_date = Sys.time(),
    tar_files_processed = length(tar_files)
)

saveRDS(extraction_summary, "docs/data_description/extraction_summary.rds")

# Create a data inventory
data_inventory <- data.frame(
    sample = extracted_samples,
    sample_type = gsub(".*_", "", extracted_samples),  # Extract sample type (NGB, IMP, IMR, etc.)
    files_count = sapply(extracted_samples, function(s) {
        length(list.files(paste0("data/raw/extracted/", s), recursive = TRUE))
    }),
    extraction_status = "success"
)

write.csv(data_inventory, "docs/data_description/data_inventory.csv", row.names = FALSE)

# Save script information
save_script_info(
    script_name = "extract_and_organize_data.R",
    description = "Extract tar.gz files and organize processed scRNA-seq data",
    output_files = list.files("data", recursive = TRUE, full.names = TRUE)
)

cat("\n=== EXTRACTION AND ORGANIZATION COMPLETE ===\n")
cat("Check docs/data_description/ for detailed information\n")
cat("Next step: Load and combine expression matrices\n")
