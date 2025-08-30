### Download Clinical Data for Gender Analysis - Simple Version
# This script downloads clinical data with gender information for GBM, OV, and PAAD

# Load required libraries
library(TCGAbiolinks)
library(dplyr)

# Function to download clinical data for a specific cancer type
download_clinical_data <- function(cancer_type, output_dir = "data/raw/clinical") {
  
  cat("=== Downloading Clinical Data for", cancer_type, "===\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define project ID
  project_id <- paste0("TCGA-", cancer_type)
  
  # Download clinical data
  cat("Downloading clinical data for", project_id, "...\n")
  
  tryCatch({
    # Get clinical data
    clinical_data <- GDCquery_clinic(project = project_id, type = "clinical")
    
    # Check if gender information is available
    gender_cols <- grep("gender|sex", colnames(clinical_data), ignore.case = TRUE, value = TRUE)
    
    if (length(gender_cols) > 0) {
      cat("✓ Gender information found in columns:", paste(gender_cols, collapse = ", "), "\n")
      
      # Show gender distribution
      for (col in gender_cols) {
        cat("Gender distribution in", col, ":\n")
        print(table(clinical_data[[col]], useNA = "ifany"))
      }
    } else {
      cat("⚠ No gender columns found in clinical data\n")
    }
    
    # Save clinical data
    output_file <- file.path(output_dir, paste0(tolower(cancer_type), "_clinical_data.rds"))
    saveRDS(clinical_data, output_file)
    
    # Also save as CSV for easy viewing
    csv_file <- file.path(output_dir, paste0(tolower(cancer_type), "_clinical_data.csv"))
    write.csv(clinical_data, csv_file, row.names = FALSE)
    
    cat("✓ Clinical data saved to:", output_file, "\n")
    cat("✓ CSV version saved to:", csv_file, "\n")
    
    return(clinical_data)
    
  }, error = function(e) {
    cat("✗ Error downloading clinical data for", cancer_type, ":", e$message, "\n")
    return(NULL)
  })
}

# Function to create gender analysis summary
create_gender_summary <- function() {
  
  cat("\n=== Creating Gender Analysis Summary ===\n")
  
  clinical_dir <- "data/raw/clinical"
  
  if (!dir.exists(clinical_dir)) {
    cat("Clinical data directory not found. Please run download_clinical_data() first.\n")
    return(NULL)
  }
  
  # Check for clinical data files
  clinical_files <- list.files(clinical_dir, pattern = "_clinical_data.rds", full.names = TRUE)
  
  if (length(clinical_files) == 0) {
    cat("No clinical data files found. Please download clinical data first.\n")
    return(NULL)
  }
  
  summary_data <- data.frame()
  
  for (file in clinical_files) {
    cancer_type <- gsub(".*/([a-z]+)_clinical_data.rds", "\\1", file)
    cancer_type <- toupper(cancer_type)
    
    cat("Processing", cancer_type, "clinical data...\n")
    
    clinical_data <- readRDS(file)
    
    # Find gender columns
    gender_cols <- grep("gender|sex", colnames(clinical_data), ignore.case = TRUE, value = TRUE)
    
    if (length(gender_cols) > 0) {
      for (col in gender_cols) {
        gender_dist <- table(clinical_data[[col]], useNA = "ifany")
        
        summary_row <- data.frame(
          Cancer_Type = cancer_type,
          Gender_Column = col,
          Male_Count = ifelse("male" %in% tolower(names(gender_dist)), 
                             gender_dist[tolower(names(gender_dist)) == "male"], 0),
          Female_Count = ifelse("female" %in% tolower(names(gender_dist)), 
                               gender_dist[tolower(names(gender_dist)) == "female"], 0),
          Unknown_Count = ifelse("unknown" %in% tolower(names(gender_dist)) | 
                                is.na(names(gender_dist)), 
                                sum(gender_dist[tolower(names(gender_dist)) == "unknown" | is.na(names(gender_dist))]), 0),
          Total_Samples = nrow(clinical_data)
        )
        
        summary_data <- rbind(summary_data, summary_row)
      }
    } else {
      summary_row <- data.frame(
        Cancer_Type = cancer_type,
        Gender_Column = "None_Found",
        Male_Count = 0,
        Female_Count = 0,
        Unknown_Count = nrow(clinical_data),
        Total_Samples = nrow(clinical_data)
      )
      summary_data <- rbind(summary_data, summary_row)
    }
  }
  
  # Save summary
  summary_file <- file.path(clinical_dir, "gender_analysis_summary.csv")
  write.csv(summary_data, summary_file, row.names = FALSE)
  
  cat("✓ Gender analysis summary saved to:", summary_file, "\n")
  
  return(summary_data)
}

# Main execution
cat("=== Clinical Data Download for Gender Analysis ===\n\n")

# Download clinical data for each cancer type
cancer_types <- c("GBM", "OV", "PAAD")

for (cancer_type in cancer_types) {
  download_clinical_data(cancer_type)
  cat("\n")
}

# Create summary
summary <- create_gender_summary()

if (!is.null(summary)) {
  cat("\n=== GENDER ANALYSIS SUMMARY ===\n")
  print(summary)
}
