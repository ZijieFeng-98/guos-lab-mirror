# TCGA Data Download Script
# Use this script when GDC server is available to download data for offline use

library(TCGAbiolinks)
library(SummarizedExperiment)

cat("=== TCGA-GBM Data Download Script ===\n\n")

# Function to check GDC server status
check_gdc_status <- function() {
  cat("Checking GDC server status...\n")
  tryCatch({
    response <- httr::GET("https://api.gdc.cancer.gov/status", timeout(10))
    if (httr::status_code(response) == 200) {
      cat("✓ GDC server is accessible\n")
      return(TRUE)
    } else {
      cat("✗ GDC server returned status code:", httr::status_code(response), "\n")
      return(FALSE)
    }
  }, error = function(e) {
    cat("✗ GDC server is not accessible:", e$message, "\n")
    return(FALSE)
  })
}

# Function to download and save data
download_and_save_data <- function() {
  cat("Querying TCGA-GBM data...\n")
  
  # Query for TCGA-GBM data
  query <- GDCquery(
    project = "TCGA-GBM",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  cat("Downloading data (this may take several minutes)...\n")
  GDCdownload(query)
  
  cat("Preparing data...\n")
  gbm_data <- GDCprepare(query)
  
  # Save data locally
  output_file <- "gbm_data.rds"
  cat("Saving data to:", output_file, "\n")
  saveRDS(gbm_data, file = output_file)
  
  cat("✓ Data successfully downloaded and saved!\n")
  cat("File size:", format(file.size(output_file), units = "MB"), "\n")
  cat("Samples:", ncol(gbm_data), "\n")
  cat("Genes:", nrow(gbm_data), "\n")
  
  return(gbm_data)
}

# Main execution
if (check_gdc_status()) {
  cat("\nProceeding with data download...\n")
  data <- download_and_save_data()
  
  cat("\n=== Download Complete ===\n")
  cat("You can now run the offline analysis script:\n")
  cat("source('TCGA_GBM_Correlation_Analysis_Offline.R')\n")
} else {
  cat("\n=== GDC Server Unavailable ===\n")
  cat("Please try again later when the server is available.\n")
  cat("You can check server status at: https://portal.gdc.cancer.gov/\n")
}
