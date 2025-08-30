# Check gene names format in OV data
library(SummarizedExperiment)

cat("=== Checking OV Gene Names Format ===\n")

# Load OV data
ov_data <- readRDS("ov_data.rds")

# Check different ways to get gene names
cat("Row names (first 10):\n")
print(head(rownames(ov_data), 10))

cat("\nRow data columns:\n")
print(colnames(rowData(ov_data)))

if ("gene_name" %in% colnames(rowData(ov_data))) {
  cat("\nGene names from rowData (first 10):\n")
  print(head(rowData(ov_data)$gene_name, 10))
}

if ("gene_id" %in% colnames(rowData(ov_data))) {
  cat("\nGene IDs from rowData (first 10):\n")
  print(head(rowData(ov_data)$gene_id, 10))
}

# Check if any of our genes exist in different formats
pathway_genes <- c("FLVCR1", "FLVCR2", "TXNDC16", "SOAT1", "SOAT2", "SCAP", "SREBF1",
                   "TMEM173", "MB21D1", "LRP1", "LRP8", "LDLR", "APOA1", "APOB", "APOC1",
                   "APOD", "APOE", "DGAT1", "DGAT2")

immune_markers <- c("IFNA1", "IFNB1", "IL1B", "IL6", "TNF", "GZMA", "GZMB", "PRF1", "CXCL9", "CXCL10")

cat("\nChecking for pathway genes in row names:\n")
for (gene in pathway_genes[1:5]) {
  found <- grep(gene, rownames(ov_data), ignore.case = TRUE)
  if (length(found) > 0) {
    cat("  ", gene, "found at positions:", found[1:3], "\n")
  } else {
    cat("  ", gene, "not found\n")
  }
}

cat("\nChecking for immune markers in row names:\n")
for (gene in immune_markers[1:5]) {
  found <- grep(gene, rownames(ov_data), ignore.case = TRUE)
  if (length(found) > 0) {
    cat("  ", gene, "found at positions:", found[1:3], "\n")
  } else {
    cat("  ", gene, "not found\n")
  }
}

