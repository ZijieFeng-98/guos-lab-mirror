#!/usr/bin/env Rscript
# qPCR Cleaner and Analyzer - Package Installation Script
# This script installs all required packages for the application

cat("Installing required packages for qPCR Cleaner and Analyzer...\n")

# List of required packages
required_packages <- c(
  "shiny",        # Web application framework
  "gtools",       # Permutation calculations
  "readxl",       # Excel file reading
  "DT",           # Interactive data tables
  "ggplot2",      # Data visualization
  "dplyr",        # Data manipulation
  "tidyr",        # Data tidying
  "shinyjs",      # JavaScript functionality for Shiny
  "shinyWidgets"  # Enhanced UI widgets
)

# Function to install packages if not already installed
install_if_missing <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", package_name, "...\n"))
    install.packages(package_name, dependencies = TRUE)
    
    # Check if installation was successful
    if (require(package_name, character.only = TRUE, quietly = TRUE)) {
      cat(paste("✓", package_name, "installed successfully\n"))
    } else {
      cat(paste("✗ Failed to install", package_name, "\n"))
      return(FALSE)
    }
  } else {
    cat(paste("✓", package_name, "already installed\n"))
  }
  return(TRUE)
}

# Install all packages
cat("\nChecking and installing packages...\n")
cat("=====================================\n")

success_count <- 0
total_packages <- length(required_packages)

for (package in required_packages) {
  if (install_if_missing(package)) {
    success_count <- success_count + 1
  }
}

cat("\n=====================================\n")
cat(paste("Installation complete:", success_count, "of", total_packages, "packages installed successfully\n"))

if (success_count == total_packages) {
  cat("\n🎉 All packages installed successfully!\n")
  cat("You can now run the qPCR Cleaner and Analyzer app with:\n")
  cat("shiny::runApp()\n")
} else {
  cat("\n⚠️  Some packages failed to install. Please check the error messages above.\n")
  cat("You may need to install them manually or check your internet connection.\n")
}

# Test if the app can be loaded
cat("\nTesting app loading...\n")
tryCatch({
  source("app.R")
  cat("✓ App file loaded successfully\n")
}, error = function(e) {
  cat("✗ Error loading app file:", e$message, "\n")
  cat("Make sure you're in the correct directory with the app.R file.\n")
})
