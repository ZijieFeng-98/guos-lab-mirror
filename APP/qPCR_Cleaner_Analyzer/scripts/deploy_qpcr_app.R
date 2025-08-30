################################################################################
# Script: deploy_qpcr_app.R
#
# This script:
# - Tests the Shiny app locally
# - Deploys the Shiny app to shinyapps.io
# - Provides deployment status and troubleshooting
#
# Author: Zijie Feng
# Last updated: 2024
################################################################################

# Load required libraries
library(shiny)
library(rsconnect)

# Configuration
APP_NAME <- "qpcr-cleaner"
APP_PATH <- getwd()  # Use current working directory

# Function to check if required packages are installed
check_packages <- function() {
  required_packages <- c("shiny", "gtools", "readxl", "DT", "ggplot2", 
                        "dplyr", "tidyr", "shinyjs", "shinyWidgets")
  
  missing_packages <- c()
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  if (length(missing_packages) > 0) {
    cat("Missing packages:", paste(missing_packages, collapse = ", "), "\n")
    cat("Installing missing packages...\n")
    install.packages(missing_packages, dependencies = TRUE)
  } else {
    cat("✓ All required packages are installed\n")
  }
}

# Function to test app locally
test_app_locally <- function() {
  cat("Testing app locally...\n")
  cat("App will open in your browser. Close the browser window to continue.\n")
  
  tryCatch({
    # Test if app.R exists
    if (!file.exists("app.R")) {
      stop("app.R file not found in current directory")
    }
    
    # Test run the app (this will open in browser)
    shiny::runApp(launch.browser = TRUE)
    cat("✓ Local test completed successfully\n")
    return(TRUE)
  }, error = function(e) {
    cat("✗ Local test failed:", e$message, "\n")
    return(FALSE)
  })
}

# Function to deploy to shinyapps.io
deploy_to_shinyapps <- function() {
  cat("Deploying to shinyapps.io...\n")
  
  tryCatch({
    # Check if rsconnect is configured
    if (!file.exists(".Rprofile")) {
      cat("Setting up rsconnect configuration...\n")
      
      # You can set your account info here or use rsconnect::setAccountInfo()
      # For security, it's better to use environment variables or .Rprofile
      cat("Please configure your shinyapps.io account:\n")
      cat("1. Go to https://www.shinyapps.io/\n")
      cat("2. Create an account and get your token\n")
      cat("3. Run: rsconnect::setAccountInfo(name='yourname', token='yourtoken', secret='yoursecret')\n")
      
      # Configure ShinyApps.io credentials
      # For security, use environment variables or .Rprofile
      cat("Please configure your ShinyApps.io credentials:\n")
      cat("1. Create a .Rprofile file in your home directory\n")
      cat("2. Add your credentials there:\n")
      cat("   rsconnect::setAccountInfo(name='your-username', token='your-token', secret='your-secret')\n")
      cat("3. Or set environment variables: SHINYAPPS_NAME, SHINYAPPS_TOKEN, SHINYAPPS_SECRET\n")
      
      return(FALSE)
    }
    
    # Deploy the app
    rsconnect::deployApp(
      appDir = APP_PATH,
      appName = APP_NAME,
      appTitle = "qPCR Cleaner and Analyzer"
    )
    
    cat("✓ Deployment successful!\n")
    cat("Your app is available at: https://zijiefeng.shinyapps.io/", APP_NAME, "\n", sep = "")
    return(TRUE)
    
  }, error = function(e) {
    cat("✗ Deployment failed:", e$message, "\n")
    cat("Troubleshooting tips:\n")
    cat("1. Check your internet connection\n")
    cat("2. Verify your shinyapps.io credentials\n")
    cat("3. Ensure all required packages are installed\n")
    cat("4. Check that app.R file exists and is valid\n")
    return(FALSE)
  })
}

# Main execution
main <- function() {
  cat("========================================\n")
  cat("qPCR Cleaner and Analyzer - Deployment\n")
  cat("========================================\n\n")
  
  # Check packages
  cat("1. Checking required packages...\n")
  check_packages()
  cat("\n")
  
  # Test locally
  cat("2. Testing app locally...\n")
  local_test <- test_app_locally()
  cat("\n")
  
  if (!local_test) {
    cat("Local test failed. Please fix any issues before deploying.\n")
    return()
  }
  
  # Deploy
  cat("3. Deploying to shinyapps.io...\n")
  deploy_success <- deploy_to_shinyapps()
  
  if (deploy_success) {
    cat("\n🎉 Deployment completed successfully!\n")
    cat("Your qPCR Cleaner and Analyzer is now live online!\n")
  } else {
    cat("\n⚠️  Deployment failed. Please check the error messages above.\n")
  }
}

# Run the main function
main()
