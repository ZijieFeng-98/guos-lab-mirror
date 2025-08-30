# Download TCGA-OV data (simple version)
library(TCGAbiolinks)
library(SummarizedExperiment)

cat("=== Downloading TCGA-OV Data ===\n")

# Step 1: Query OV data
cat("Querying TCGA-OV data...\n")
tryCatch({
  query <- GDCquery(
    project = "TCGA-OV",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  cat("Found", length(query$results), "files to download\n")
  
  # Step 2: Download data
  cat("Downloading data...\n")
  GDCdownload(query)
  
  # Step 3: Prepare data
  cat("Preparing data...\n")
  ov_data <- GDCprepare(query)
  
  cat("✓ Data prepared successfully\n")
  cat("  Genes:", nrow(ov_data), "\n")
  cat("  Samples:", ncol(ov_data), "\n")
  
  # Step 4: Save as RDS
  cat("Saving to ov_data.rds...\n")
  saveRDS(ov_data, file = "ov_data.rds")
  
  cat("✓ OV data saved successfully!\n")
  cat("  File size:", format(file.size("ov_data.rds"), units = "MB"), "\n")
  cat("  Genes:", nrow(ov_data), "\n")
  cat("  Samples:", ncol(ov_data), "\n")
  
  # Step 5: Print sample information
  cat("\nSample information:\n")
  sample_types <- table(colData(ov_data)$shortLetterCode)
  for (type in names(sample_types)) {
    cat("  ", type, ":", sample_types[type], "samples\n")
  }
  
  cat("\n✓ Download and conversion complete!\n")
  cat("You can now run: TCGA_OV_Correlation_Analysis_Offline.R\n")
  
}, error = function(e) {
  cat("✗ Error:", e$message, "\n")
  cat("\n=== GDC Server Issues - Manual Download Instructions ===\n")
  cat("1. Visit: https://portal.gdc.cancer.gov/\n")
  cat("2. Navigate to 'Repository' -> 'Files'\n")
  cat("3. Apply filters:\n")
  cat("   - Project: TCGA-OV\n")
  cat("   - Data Category: Transcriptome Profiling\n")
  cat("   - Data Type: Gene Expression Quantification\n")
  cat("   - Workflow Type: STAR - Counts\n")
  cat("4. Download the manifest file\n")
  cat("5. Use GDC Data Transfer Tool or gdc-client\n")
  cat("6. Convert to RDS format using the conversion scripts\n")
})

