#!/usr/bin/env Rscript

# ==============================================================================
# GSE222520 Analysis Template
# Template for single-cell RNA-seq analysis
# ==============================================================================

# Load required libraries
if (!requireNamespace("Seurat", quietly = TRUE))
    install.packages("Seurat")

library(Seurat)
library(dplyr)
library(ggplot2)

# Load the combined data
cat("Loading combined expression matrix...\n")
expression_matrix <- readRDS("data/processed/GSE222520_combined_expression_matrix.rds")
cell_metadata <- read.csv("data/processed/GSE222520_cell_metadata.csv")
gene_metadata <- read.csv("data/processed/GSE222520_gene_metadata.csv")

# Create Seurat object
cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(counts = expression_matrix, meta.data = cell_metadata)

# Add sample information
seurat_obj$sample_type <- cell_metadata$sample_type
seurat_obj$sample <- cell_metadata$sample

cat("Seurat object created with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n")

# Basic QC
cat("\n=== QUALITY CONTROL ===\n")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Plot QC metrics
pdf("figures/qc/qc_metrics.pdf", width = 12, height = 8)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

cat("QC plots saved to figures/qc/qc_metrics.pdf\n")

# Save Seurat object
saveRDS(seurat_obj, "data/processed/GSE222520_seurat_object.rds")
cat("Seurat object saved to data/processed/GSE222520_seurat_object.rds\n")

cat("\n=== ANALYSIS TEMPLATE COMPLETE ===\n")
cat("Ready for further analysis steps\n")
