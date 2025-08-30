# TCGA OV Correlation Analysis: Expanded Immune Marker List (Offline Version)
# High-grade Serous Ovarian Cancer Analysis using local RDS file

# Load required libraries
library(SummarizedExperiment)
library(pheatmap)

cat("=== TCGA OV Correlation Analysis (Offline Mode) ===\n")

# Function to load data from local file
load_local_data <- function(file_path) {
  if (file.exists(file_path)) {
    cat("Loading data from local file:", file_path, "\n")
    return(readRDS(file_path))
  } else {
    stop("Local data file not found. Please convert OV data to RDS format first.")
  }
}

# Main analysis function
run_ov_analysis <- function(data_file = "ov_data.rds") {
  # Load data
  ov_data <- load_local_data(data_file)
  
  cat("✓ Data loaded successfully\n")
  cat("  Genes:", nrow(ov_data), "\n")
  cat("  Samples:", ncol(ov_data), "\n")
  
  # ✅ Keep only primary tumor samples (shortLetterCode == 'TP')
  if ("shortLetterCode" %in% colnames(colData(ov_data))) {
    ov_data <- ov_data[, colData(ov_data)$shortLetterCode == "TP"]
    cat("  Primary tumor samples:", ncol(ov_data), "\n")
  } else {
    cat("  Note: No shortLetterCode column found, using all samples\n")
  }
  
  # Step 2: Extract expression matrix and clean gene names
  expr <- assay(ov_data)
  gene_names <- rowData(ov_data)$gene_name
  valid_genes <- !is.na(gene_names) & !duplicated(gene_names)
  expr <- expr[valid_genes, ]
  rownames(expr) <- gene_names[valid_genes]
  
  # Step 3: Define genes of interest
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
  
  # Step 4: Filter expression matrix
  genes_of_interest <- intersect(rownames(expr), unique(c(pathway_genes, immune_markers)))
  expr_filtered <- expr[genes_of_interest, ]
  
  # Step 5: Spearman correlation
  cor_mat <- cor(t(expr_filtered), method = "spearman")
  pathway_in_expr <- intersect(pathway_genes, rownames(expr_filtered))
  immune_in_expr <- intersect(immune_markers, rownames(expr_filtered))
  cor_block <- cor_mat[pathway_in_expr, immune_in_expr]
  
  # Step 6: Create output directory if it doesn't exist
  out_dir <- file.path("Processed_Data")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Step 7: Save heatmap
  out_path <- file.path(out_dir, "Spearman_Correlation_Heatmap_OV_All_Immune.png")
  
  png(filename = out_path, width = 2600, height = 1600, res = 300)
  pheatmap(cor_block,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = FALSE,
           fontsize = 10,
           main = "Spearman Correlation: Pathway vs Full Immune Marker Set (TCGA-OV)")
  dev.off()
  
  # Step 8: Print summary statistics
  cat("Analysis Summary (TCGA-OV):\n")
  cat("===========================\n")
  cat("Total samples analyzed:", ncol(expr_filtered), "\n")
  cat("Pathway genes found in data:", length(pathway_in_expr), "/", length(pathway_genes), "\n")
  cat("Immune markers found in data:", length(immune_in_expr), "/", length(immune_markers), "\n")
  cat("Output saved to:", out_path, "\n")
  
  # Step 9: Save correlation matrix as CSV for further analysis
  cor_csv_path <- file.path(out_dir, "Spearman_Correlation_Matrix_OV_All_Immune.csv")
  write.csv(cor_block, file = cor_csv_path)
  
  cat("Correlation matrix saved to:", cor_csv_path, "\n")
  
  return(list(correlation_matrix = cor_block, 
              pathway_genes = pathway_in_expr, 
              immune_genes = immune_in_expr))
}

# Run the analysis
tryCatch({
  result <- run_ov_analysis()
  cat("\n✓ Analysis completed successfully!\n")
}, error = function(e) {
  cat("\n✗ Error:", e$message, "\n")
  cat("\n=== OV Data Not Available - Manual Instructions ===\n")
  cat("The OV data files appear to be corrupted. You have two options:\n\n")
  cat("1. Re-download OV data from GDC:\n")
  cat("   - Visit: https://portal.gdc.cancer.gov/\n")
  cat("   - Navigate to 'Repository' -> 'Files'\n")
  cat("   - Apply filters: Project=TCGA-OV, Data Category=Transcriptome Profiling\n")
  cat("   - Download and convert to RDS format\n\n")
  cat("2. Continue with available data:\n")
  cat("   - You have GBM and BRCA data available\n")
  cat("   - Run: TCGA_GBM_Correlation_Analysis_Offline.R\n")
  cat("   - Run: TCGA_BRCA_Correlation_Analysis_Offline.R\n")
})
