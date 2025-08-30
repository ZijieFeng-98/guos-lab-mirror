# Package Installation Script for TCGA GBM Correlation Analysis
# Run this script first to ensure all required packages are installed

cat("Installing and checking required packages for TCGA GBM Correlation Analysis...\n\n")

# List of required packages
required_packages <- c("TCGAbiolinks", "SummarizedExperiment", "pheatmap")

# Function to install package if not already installed
install_if_missing <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", package_name, "...\n")
    install.packages(package_name, dependencies = TRUE)
    if (require(package_name, character.only = TRUE, quietly = TRUE)) {
      cat("✓", package_name, "installed successfully\n")
    } else {
      cat("✗ Failed to install", package_name, "\n")
      return(FALSE)
    }
  } else {
    cat("✓", package_name, "already installed\n")
  }
  return(TRUE)
}

# Install all required packages
cat("Checking and installing required packages:\n")
cat("==========================================\n")

all_installed <- TRUE
for (pkg in required_packages) {
  if (!install_if_missing(pkg)) {
    all_installed <- FALSE
  }
}

cat("\n")
if (all_installed) {
  cat("✓ All required packages are ready!\n")
  cat("You can now run: source('TCGA_GBM_Correlation_Analysis_Expanded_Immune.R')\n")
} else {
  cat("✗ Some packages failed to install. Please check your R installation and internet connection.\n")
}

# Check R version
cat("\nR version information:\n")
cat("=====================\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")

# Memory check
cat("\nMemory information:\n")
cat("==================\n")
if (require(pryr, quietly = TRUE)) {
  cat("Available memory:", format(pryr::mem_used(), units = "MB"), "\n")
} else {
  cat("Install 'pryr' package for detailed memory information\n")
}
