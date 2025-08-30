# Mouse Management Software - Package Installation Script
# Run this script to install all required packages

required_packages <- c(
  "shiny",           # Web application framework
  "shinydashboard",  # Dashboard UI components
  "DT",              # DataTables for interactive tables
  "dplyr",           # Data manipulation
  "readxl",          # Read Excel files
  "writexl",         # Write Excel files
  "shinyjs",         # JavaScript functionality
  "htmltools",       # HTML generation
  "jsonlite"         # JSON handling
)

# Install missing packages
cat("Installing required packages...\n")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(pkg, "✓\n")
  }
}

# Verify installation
cat("\nVerifying installation...\n")
missing <- c()
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing <- c(missing, pkg)
  }
}

if (length(missing) == 0) {
  cat("✅ All packages installed successfully!\n")
  cat("Run: source('app.R') to start the application\n")
} else {
  cat("❌ Failed to install:", paste(missing, collapse = ", "), "\n")
}
