#!/usr/bin/env Rscript
# qPCR Cleaner and Analyzer - Launcher Script
# Simple script to run the app

cat("Starting qPCR Cleaner and Analyzer...\n")

# Check if required packages are installed
required_packages <- c("shiny", "gtools", "readxl", "DT", "ggplot2", 
                      "dplyr", "tidyr", "shinyjs", "shinyWidgets")

missing_packages <- c()
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cran.rstudio.com/")
}

# Run the app
cat("Launching app in browser...\n")
shiny::runApp(launch.browser = TRUE)
