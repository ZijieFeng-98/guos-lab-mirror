### TCGA GBM Correlation Analysis: Expanded Immune Marker List (from Merged_Immune_Pathway_Markers_in_Glioma.docx)

# Load required libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(pheatmap)

# Step 1: Query and download data
query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
gbm_data <- GDCprepare(query)

# ✅ Keep only primary tumor samples (shortLetterCode == 'TP')
gbm_data <- gbm_data[, colData(gbm_data)$shortLetterCode == "TP"]

# Step 2: Extract expression matrix and clean gene names
expr <- assay(gbm_data)
gene_names <- rowData(gbm_data)$gene_name
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

# Step 6: Save heatmap
out_path <- file.path("/Users/zijiefeng/Desktop/Guo's lab",
                      "My_Research",
                      "Immunology_analysis",
                      "Processed_Data",
                      "Spearman_Correlation_Heatmap_All_Immune.png")

png(filename = out_path, width = 2600, height = 1600, res = 300)
pheatmap(cor_block,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = FALSE,
         fontsize = 10,
         main = "Spearman Correlation: Pathway vs Full Immune Marker Set (TCGA-GBM)")
dev.off()
