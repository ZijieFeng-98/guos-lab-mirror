# GDC Server Status Checker
# Quick script to check if GDC server is available

cat("=== GDC Server Status Check ===\n\n")

# Check if httr package is available
if (!require(httr, quietly = TRUE)) {
  cat("Installing httr package...\n")
  install.packages("httr")
  library(httr)
}

# Function to check GDC status
check_gdc_status <- function() {
  cat("Checking GDC server status...\n")
  
  tryCatch({
    # Set timeout to 10 seconds
    response <- httr::GET("https://api.gdc.cancer.gov/status", 
                         httr::timeout(10))
    
    if (httr::status_code(response) == 200) {
      cat("✓ GDC server is accessible!\n")
      cat("Status code:", httr::status_code(response), "\n")
      return(TRUE)
    } else {
      cat("✗ GDC server returned error status:", httr::status_code(response), "\n")
      return(FALSE)
    }
  }, error = function(e) {
    cat("✗ GDC server is not accessible\n")
    cat("Error:", e$message, "\n")
    return(FALSE)
  })
}

# Check status
server_available <- check_gdc_status()

cat("\n=== Recommendations ===\n")

if (server_available) {
  cat("✓ Server is available! You can:\n")
  cat("1. Run the main analysis: source('TCGA_GBM_Correlation_Analysis_Expanded_Immune.R')\n")
  cat("2. Download data for offline use: source('download_tcga_data.R')\n")
} else {
  cat("✗ Server is down. You can:\n")
  cat("1. Check if you have local data: ls -la gbm_data.rds\n")
  cat("2. Run offline analysis if data exists: source('TCGA_GBM_Correlation_Analysis_Offline.R')\n")
  cat("3. Wait and try again later\n")
  cat("4. Visit https://portal.gdc.cancer.gov/ for manual download\n")
}

cat("\n=== Current Directory Files ===\n")
cat("Available scripts:\n")
files <- list.files(pattern = "\\.R$")
for (file in files) {
  cat("-", file, "\n")
}

if (file.exists("gbm_data.rds")) {
  cat("\n✓ Local data file found: gbm_data.rds\n")
  file_info <- file.info("gbm_data.rds")
  cat("File size:", format(file_info$size, units = "MB"), "\n")
  cat("Last modified:", file_info$mtime, "\n")
} else {
  cat("\n✗ No local data file found\n")
}
