### TCGA PAAD Correlation Analysis: Offline Version (Using Downloaded Data)
# Pancreatic Adenocarcinoma Analysis

# Load required libraries
library(SummarizedExperiment)
library(pheatmap)

# Function to load data from local file
load_local_data <- function(file_path) {
  if (file.exists(file_path)) {
    cat("Loading data from local file:", file_path, "\n")
    return(readRDS(file_path))
  } else {
    stop("Local data file not found. Please download data manually or wait for GDC server to be available.")
  }
}

# Main analysis function
run_paad_offline_analysis <- function(data_file = "paad_data.rds") {
  
  cat("=== TCGA PAAD Correlation Analysis (Offline Mode) ===\n\n")
  
  # Check if we have local data
  if (file.exists(data_file)) {
    cat("✓ Found local PAAD data file\n")
    paad_data <- load_local_data(data_file)
  } else {
    cat("✗ No local PAAD data found\n")
    cat("Please download TCGA-PAAD data and save as 'paad_data.rds'\n")
    return(NULL)
  }
  
  # Continue with analysis if data is available
  cat("Proceeding with PAAD analysis...\n\n")
  
  # Step 1: Filter for primary tumor samples
  if ("shortLetterCode" %in% colnames(colData(paad_data))) {
    paad_data <- paad_data[, colData(paad_data)$shortLetterCode == "TP"]
    cat("✓ Filtered for primary tumor samples:", ncol(paad_data), "samples\n")
  } else {
    cat("⚠ Warning: No sample type information found\n")
  }
  
  # Step 2: Extract expression matrix and clean gene names
  expr <- assay(paad_data)
  gene_names <- rowData(paad_data)$gene_name
  
  if (is.null(gene_names)) {
    # Try alternative gene name columns
    if ("gene_name" %in% colnames(rowData(paad_data))) {
      gene_names <- rowData(paad_data)$gene_name
    } else if ("gene_id" %in% colnames(rowData(paad_data))) {
      gene_names <- rowData(paad_data)$gene_id
    } else {
      gene_names <- rownames(expr)
    }
  }
  
  valid_genes <- !is.na(gene_names) & !duplicated(gene_names)
  expr <- expr[valid_genes, ]
  rownames(expr) <- gene_names[valid_genes]
  
  cat("✓ Expression matrix prepared:", nrow(expr), "genes,", ncol(expr), "samples\n")
  
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
  
  cat("✓ Genes of interest found:", length(genes_of_interest), "\n")
  cat("  - Pathway genes:", length(intersect(pathway_genes, genes_of_interest)), "/", length(pathway_genes), "\n")
  cat("  - Immune markers:", length(intersect(immune_markers, genes_of_interest)), "/", length(immune_markers), "\n")
  
  # Step 5: Spearman correlation
  cat("\nComputing correlations...\n")
  cor_mat <- cor(t(expr_filtered), method = "spearman")
  pathway_in_expr <- intersect(pathway_genes, rownames(expr_filtered))
  immune_in_expr <- intersect(immune_markers, rownames(expr_filtered))
  cor_block <- cor_mat[pathway_in_expr, immune_in_expr]
  
  # Step 6: Create output directory
  out_dir <- file.path("/Users/zijiefeng/Desktop/Guo's lab",
                       "My_Research",
                       "Immunology_analysis",
                       "Processed_Data")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Step 7: Save heatmap
  out_path <- file.path(out_dir, "Spearman_Correlation_Heatmap_PAAD_All_Immune_Offline.png")
  
  png(filename = out_path, width = 2600, height = 1600, res = 300)
  pheatmap(cor_block,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = FALSE,
           fontsize = 10,
           main = "Spearman Correlation: Pathway vs Full Immune Marker Set (TCGA-PAAD) - Offline Analysis")
  dev.off()
  
  # Step 8: Save correlation matrix
  cor_csv_path <- file.path(out_dir, "Spearman_Correlation_Matrix_PAAD_All_Immune_Offline.csv")
  write.csv(cor_block, file = cor_csv_path)
  
  # Step 9: Print summary
  cat("\n=== PAAD Analysis Complete ===\n")
  cat("Total samples analyzed:", ncol(expr_filtered), "\n")
  cat("Pathway genes found:", length(pathway_in_expr), "/", length(pathway_genes), "\n")
  cat("Immune markers found:", length(immune_in_expr), "/", length(immune_markers), "\n")
  cat("Output saved to:", out_path, "\n")
  cat("Correlation matrix saved to:", cor_csv_path, "\n")
  
  return(list(correlation_matrix = cor_block, 
              pathway_genes = pathway_in_expr, 
              immune_genes = immune_in_expr))
}

# Run the offline analysis
result <- run_paad_offline_analysis()
