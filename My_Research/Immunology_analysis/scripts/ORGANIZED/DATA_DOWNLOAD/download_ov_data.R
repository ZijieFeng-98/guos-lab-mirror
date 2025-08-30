# Download TCGA-OV data and convert to RDS format
library(TCGAbiolinks)
library(SummarizedExperiment)

cat("=== Downloading TCGA-OV Data ===\n")

# Check GDC server status
cat("Checking GDC server status...\n")
tryCatch({
  status <- GDCstatus()
  cat("✓ GDC server is accessible\n")
}, error = function(e) {
  cat("✗ GDC server error:", e$message, "\n")
  cat("Please try again later or use manual download\n")
  stop("GDC server unavailable")
})

# Step 1: Query OV data
cat("Querying TCGA-OV data...\n")
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

