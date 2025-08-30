# Check genes in OV data
library(SummarizedExperiment)

cat("=== Checking OV Data Genes ===\n")

# Load OV data
ov_data <- readRDS("ov_data.rds")
cat("✓ OV data loaded\n")
cat("  Genes:", nrow(ov_data), "\n")
cat("  Samples:", ncol(ov_data), "\n")

# Filter for primary tumor samples
ov_data <- ov_data[, colData(ov_data)$shortLetterCode == "TP"]
cat("  Primary tumor samples:", ncol(ov_data), "\n")

# Extract expression matrix
expr <- assay(ov_data)
gene_names <- rownames(ov_data)
valid_genes <- !is.na(gene_names) & !duplicated(gene_names)
expr <- expr[valid_genes, ]
rownames(expr) <- gene_names[valid_genes]

cat("  Valid genes:", nrow(expr), "\n")

# Define genes of interest
pathway_genes <- c("FLVCR1", "FLVCR2", "TXNDC16", "SOAT1", "SOAT2", "SCAP", "SREBF1",
                   "TMEM173", "MB21D1", "LRP1", "LRP8", "LDLR", "APOA1", "APOB", "APOC1",
                   "APOD", "APOE", "DGAT1", "DGAT2")

immune_markers <- unique(c(
  # Cytokines / Cytolytic
  "IFNA1", "IFNB1", "IL1B", "IL6", "TNF", "GZMA", "GZMB", "PRF1", "CXCL9", "CXCL10", "CCL2", "CCL4", "CCL5",
  # IFN-related
  "STAT1", "IRF7", "ISG15", "IFIT1", "IFI6", "MX1", "OAS1",
  # Immune checkpoints
  "CD274", "HAVCR2", "IDO1", "LGALS9",
  # MHC / Antigen Presentation
  "HLA-DRA", "HLA-DRB1", "CD74", "CD83", "CTSB", "CTSS",
  # T cell markers
  "CD3D", "CD4", "CD8A",
  # B cell
  "CD19", "CD79A",
  # Macrophage / Monocyte
  "CD14", "CD68", "CD163", "MRC1", "ITGAM", "CCR2", "FCN1", "LYZ",
  # Microglia
  "P2RY12", "TMEM119", "CX3CR1", "TREM2", "GPR34",
  # DC markers
  "CD1C", "ITGAX", "CLEC9A", "CD86",
  # NK cell markers
  "NCAM1", "NKG7", "KLRD1", "GNLY"
))

# Check which genes are found
pathway_in_expr <- intersect(pathway_genes, rownames(expr))
immune_in_expr <- intersect(immune_markers, rownames(expr))

cat("\nGene Analysis:\n")
cat("==============\n")
cat("Pathway genes found:", length(pathway_in_expr), "/", length(pathway_genes), "\n")
cat("Immune markers found:", length(immune_in_expr), "/", length(immune_markers), "\n")

if (length(pathway_in_expr) > 0) {
  cat("\nPathway genes found:\n")
  cat(paste(pathway_in_expr, collapse = ", "), "\n")
}

if (length(immune_in_expr) > 0) {
  cat("\nImmune markers found (first 10):\n")
  cat(paste(head(immune_in_expr, 10), collapse = ", "), "\n")
}

# Check if we have enough genes for correlation
total_genes <- length(pathway_in_expr) + length(immune_in_expr)
cat("\nTotal genes for correlation:", total_genes, "\n")

if (total_genes < 2) {
  cat("✗ Not enough genes for correlation analysis\n")
} else {
  cat("✓ Sufficient genes for correlation analysis\n")
}

