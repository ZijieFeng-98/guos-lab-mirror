#!/usr/bin/env Rscript
# GitHub Setup Script for qPCR Cleaner and Analyzer
# This script helps set up the repository for GitHub

cat("Setting up qPCR Cleaner and Analyzer for GitHub...\n")

# Function to update placeholder text in files
update_placeholders <- function() {
  cat("Updating placeholder text in files...\n")
  
  # Get GitHub username from user
  cat("Please enter your GitHub username: ")
  github_username <- readline()
  
  if (github_username == "") {
    github_username <- "yourusername"
    cat("Using default username 'yourusername'. Please update manually.\n")
  }
  
  # Files to update
  files_to_update <- c(
    "README.md",
    "DESCRIPTION",
    "deploy_qpcr_app.R"
  )
  
  for (file in files_to_update) {
    if (file.exists(file)) {
      # Read file content
      content <- readLines(file, warn = FALSE)
      
      # Replace placeholders
      content <- gsub("yourusername", github_username, content)
      content <- gsub("your-email@example.com", paste0(github_username, "@example.com"), content)
      
      # Write back to file
      writeLines(content, file)
      cat("✓ Updated", file, "\n")
    }
  }
}

# Function to check Git setup
check_git_setup <- function() {
  cat("Checking Git setup...\n")
  
  # Check if Git is installed
  git_check <- system("git --version", intern = TRUE, ignore.stderr = TRUE)
  if (length(git_check) == 0) {
    cat("✗ Git is not installed. Please install Git first.\n")
    return(FALSE)
  } else {
    cat("✓ Git is installed:", git_check[1], "\n")
  }
  
  # Check if this is a Git repository
  if (!dir.exists(".git")) {
    cat("Initializing Git repository...\n")
    system("git init")
    cat("✓ Git repository initialized\n")
  } else {
    cat("✓ Git repository already exists\n")
  }
  
  return(TRUE)
}

# Function to create initial commit
create_initial_commit <- function() {
  cat("Creating initial commit...\n")
  
  # Add all files
  system("git add .")
  
  # Create commit
  commit_result <- system("git commit -m 'Initial commit: qPCR Cleaner and Analyzer'", 
                         intern = TRUE, ignore.stderr = TRUE)
  
  if (length(commit_result) > 0) {
    cat("✓ Initial commit created\n")
  } else {
    cat("⚠️  No changes to commit or commit failed\n")
  }
}

# Function to provide GitHub setup instructions
provide_github_instructions <- function() {
  cat("\n========================================\n")
  cat("GitHub Setup Instructions\n")
  cat("========================================\n\n")
  
  cat("1. Create a new repository on GitHub:\n")
  cat("   - Go to https://github.com/new\n")
  cat("   - Name: qPCR_Cleaner_Analyzer\n")
  cat("   - Description: qPCR Cleaner and Analyzer - A comprehensive Shiny app for qPCR data analysis\n")
  cat("   - Make it Public\n")
  cat("   - Don't initialize with README (we already have one)\n\n")
  
  cat("2. Connect your local repository to GitHub:\n")
  cat("   git remote add origin https://github.com/YOUR_USERNAME/qPCR_Cleaner_Analyzer.git\n")
  cat("   git branch -M main\n")
  cat("   git push -u origin main\n\n")
  
  cat("3. Set up GitHub Actions (optional):\n")
  cat("   - Go to your repository Settings > Secrets and variables > Actions\n")
  cat("   - Add the following secrets:\n")
  cat("     * SHINYAPPS_NAME: Your ShinyApps.io username\n")
  cat("     * SHINYAPPS_TOKEN: Your ShinyApps.io token\n")
  cat("     * SHINYAPPS_SECRET: Your ShinyApps.io secret\n\n")
  
  cat("4. Enable GitHub Pages (optional):\n")
  cat("   - Go to Settings > Pages\n")
  cat("   - Source: Deploy from a branch\n")
  cat("   - Branch: main, folder: / (root)\n\n")
  
  cat("5. Update repository description and topics:\n")
  cat("   - Add topics: r, shiny, qpcr, bioinformatics, data-analysis\n")
  cat("   - Update description if needed\n\n")
}

# Main execution
main <- function() {
  cat("========================================\n")
  cat("qPCR Cleaner and Analyzer - GitHub Setup\n")
  cat("========================================\n\n")
  
  # Update placeholders
  update_placeholders()
  cat("\n")
  
  # Check Git setup
  git_ok <- check_git_setup()
  cat("\n")
  
  if (git_ok) {
    # Create initial commit
    create_initial_commit()
    cat("\n")
  }
  
  # Provide instructions
  provide_github_instructions()
  
  cat("🎉 Setup complete! Follow the instructions above to complete GitHub setup.\n")
}

# Run the main function
main()
