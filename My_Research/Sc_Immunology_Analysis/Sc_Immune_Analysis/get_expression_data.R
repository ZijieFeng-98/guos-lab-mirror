#!/usr/bin/env Rscript

# Get actual expression data for GSE222520
library(GEOquery)
library(Biobase)

# Load the saved data
gse <- readRDS("data/raw/GSE222520_raw.rds")
eset <- gse[[1]]

# Check phenotype data for supplementary file information
pheno_data <- pData(eset)
cat("Supplementary files available:\n")
print(pheno_data$supplementary_file_1)

# Check if there are any supplementary files
supp_files <- pheno_data$supplementary_file_1
supp_files <- supp_files[!is.na(supp_files) & supp_files != ""]

if (length(supp_files) > 0) {
    cat("\nFound supplementary files. Attempting to download...\n")
    
    # Download supplementary files
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
} else {
    cat("\nNo supplementary files found in phenotype data.\n")
}

# Let's also check if we can get the data directly from GEO
cat("\nAttempting to get data directly from GEO...\n")
tryCatch({
    # Try to get the full dataset
    gse_full <- getGEO("GSE222520", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/raw")
    
    if (length(gse_full) > 0) {
        eset_full <- gse_full[[1]]
        expr_matrix <- exprs(eset_full)
        
        if (nrow(expr_matrix) > 0) {
            cat("Found expression matrix with dimensions:", dim(expr_matrix), "\n")
            saveRDS(expr_matrix, "data/processed/GSE222520_expression_matrix_full.rds")
            write.csv(expr_matrix, "data/processed/GSE222520_expression_matrix_full.csv")
            cat("Saved full expression matrix\n")
        } else {
            cat("Expression matrix is still empty\n")
        }
    }
}, error = function(e) {
    cat("Error getting full dataset:", e$message, "\n")
})

# Check what files we have now
cat("\nCurrent files in data/raw:\n")
list.files("data/raw", full.names = TRUE)

cat("\nCurrent files in data/processed:\n")
list.files("data/processed", full.names = TRUE)
