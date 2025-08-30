#!/usr/bin/env Rscript

# Single-Cell RNA-seq Data Cleaning and Organization Script
# GSE222520 - Brain Tumor Leukocyte Single-Cell Analysis
# Author: Data Analysis Pipeline
# Date: 2024

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(tidyr)
library(stringr)
library(Matrix)

# Set random seed for reproducibility
set.seed(42)

# Create output directories
dir.create("data/processed/expression_matrices", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/cell_metadata", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/gene_metadata", recursive = TRUE, showWarnings = FALSE)
dir.create("results/qc_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/cleaned_data", recursive = TRUE, showWarnings = FALSE)

# Function to load 10x data
load_10x_data <- function(data_path, sample_name) {
  cat("Loading data for sample:", sample_name, "\n")
  
  # Check if filtered_feature_bc_matrix exists
  matrix_path <- file.path(data_path, "filtered_feature_bc_matrix")
  if (!dir.exists(matrix_path)) {
    # Look for the actual matrix directory
    subdirs <- list.dirs(data_path, recursive = FALSE)
    matrix_path <- file.path(subdirs[1], "filtered_feature_bc_matrix")
  }
  
  if (!dir.exists(matrix_path)) {
    stop(paste("Matrix directory not found for sample:", sample_name))
  }
  
  # Load the data
  data <- Read10X(matrix_path)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = sample_name,
    min.cells = 3,  # Genes detected in at least 3 cells
    min.features = 200  # Cells with at least 200 genes
  )
  
  # Add sample metadata
  seurat_obj$sample_id <- sample_name
  
  return(seurat_obj)
}

# Function to extract sample information from metadata
extract_sample_info <- function(metadata_file) {
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  
  # Extract key information
  sample_info <- data.frame(
    sample_id = metadata$title,
    geo_accession = metadata$geo_accession,
    tissue = metadata$`tissue:ch1`,
    diagnosis = metadata$`diagnosis:ch1`,
    genotype = metadata$`genotype:ch1`,
    primary_recurrent = metadata$`primary/recurrent:ch1`,
    treatment = metadata$`treatment:ch1`,
    cell_type = metadata$`cell type:ch1`,
    stringsAsFactors = FALSE
  )
  
  # Clean up sample names
  sample_info$sample_id <- gsub(", replicate1, scRNAseq", "", sample_info$sample_id)
  
  return(sample_info)
}

# Load sample metadata
cat("Loading sample metadata...\n")
sample_info <- extract_sample_info("data/metadata/GSE222520_phenotype_data.csv")

# Define sample groups
sample_info$group <- case_when(
  grepl("^NGB", sample_info$sample_id) ~ "Normal",
  grepl("^IMP", sample_info$sample_id) ~ "Primary_IDH_Mutant",
  grepl("^IWP", sample_info$sample_id) ~ "Primary_IDH_WildType",
  grepl("^IMR", sample_info$sample_id) ~ "Recurrent_IDH_Mutant",
  grepl("^IWR", sample_info$sample_id) ~ "Recurrent_IDH_WildType"
)

# List of sample directories
sample_dirs <- list.dirs("data/raw/extracted", recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[!grepl("^\\.", sample_dirs)]  # Remove hidden files

cat("Found", length(sample_dirs), "sample directories\n")

# Load all samples
seurat_objects <- list()
sample_paths <- list()

for (sample_dir in sample_dirs) {
  sample_name <- gsub("GSM[0-9]+_", "", sample_dir)
  sample_path <- file.path("data/raw/extracted", sample_dir)
  
  tryCatch({
    seurat_obj <- load_10x_data(sample_path, sample_name)
    seurat_objects[[sample_name]] <- seurat_obj
    sample_paths[[sample_name]] <- sample_path
    cat("Successfully loaded:", sample_name, "with", ncol(seurat_obj), "cells\n")
  }, error = function(e) {
    cat("Error loading", sample_name, ":", e$message, "\n")
  })
}

# Quality Control and Filtering
cat("\n=== QUALITY CONTROL AND FILTERING ===\n")

# Function to perform QC on individual samples
perform_qc <- function(seurat_obj, sample_name) {
  cat("Performing QC for:", sample_name, "\n")
  
  # Calculate QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  
  # Create QC plots
  p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0.1)
  p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  # Save QC plots
  ggsave(filename = file.path("results/qc_plots", paste0(sample_name, "_qc_violin.pdf")), 
         plot = p1, width = 12, height = 4)
  ggsave(filename = file.path("results/qc_plots", paste0(sample_name, "_qc_scatter.pdf")), 
         plot = p2 + p3, width = 10, height = 5)
  
  # Print QC summary
  cat("QC Summary for", sample_name, ":\n")
  cat("  Total cells:", ncol(seurat_obj), "\n")
  cat("  Total genes:", nrow(seurat_obj), "\n")
  cat("  Median genes per cell:", median(seurat_obj$nFeature_RNA), "\n")
  cat("  Median UMI per cell:", median(seurat_obj$nCount_RNA), "\n")
  cat("  Median % MT:", median(seurat_obj$percent.mt), "\n")
  
  return(seurat_obj)
}

# Apply QC to all samples
seurat_objects_qc <- list()
for (sample_name in names(seurat_objects)) {
  seurat_objects_qc[[sample_name]] <- perform_qc(seurat_objects[[sample_name]], sample_name)
}

# Filtering criteria
cat("\n=== APPLYING FILTERS ===\n")

# Function to filter samples based on QC metrics
filter_sample <- function(seurat_obj, sample_name) {
  cat("Filtering", sample_name, "...\n")
  
  # Get initial cell count
  initial_cells <- ncol(seurat_obj)
  
  # Apply filters
  seurat_obj <- subset(seurat_obj, 
                       nFeature_RNA > 200 & 
                       nFeature_RNA < 6000 & 
                       nCount_RNA > 500 & 
                       nCount_RNA < 50000 & 
                       percent.mt < 20)
  
  # Get final cell count
  final_cells <- ncol(seurat_obj)
  
  cat("  Cells before filtering:", initial_cells, "\n")
  cat("  Cells after filtering:", final_cells, "\n")
  cat("  Cells removed:", initial_cells - final_cells, "(", 
      round((initial_cells - final_cells)/initial_cells * 100, 1), "%)\n")
  
  return(seurat_obj)
}

# Apply filtering to all samples
seurat_objects_filtered <- list()
for (sample_name in names(seurat_objects_qc)) {
  seurat_objects_filtered[[sample_name]] <- filter_sample(seurat_objects_qc[[sample_name]], sample_name)
}

# Data Integration
cat("\n=== DATA INTEGRATION ===\n")

# Merge all samples
cat("Merging all samples...\n")
combined_seurat <- merge(seurat_objects_filtered[[1]], 
                        y = seurat_objects_filtered[-1], 
                        add.cell.ids = names(seurat_objects_filtered),
                        project = "GSE222520_Brain_Tumor")

# Add sample metadata
sample_meta <- sample_info[match(combined_seurat$sample_id, sample_info$sample_id), ]
combined_seurat$tissue <- sample_meta$tissue
combined_seurat$diagnosis <- sample_meta$diagnosis
combined_seurat$genotype <- sample_meta$genotype
combined_seurat$primary_recurrent <- sample_meta$primary_recurrent
combined_seurat$treatment <- sample_meta$treatment
combined_seurat$group <- sample_meta$group

# Normalize and scale data
cat("Normalizing and scaling data...\n")
combined_seurat <- NormalizeData(combined_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 2000)
combined_seurat <- ScaleData(combined_seurat, features = rownames(combined_seurat))

# Dimensionality reduction
cat("Performing dimensionality reduction...\n")
combined_seurat <- RunPCA(combined_seurat, features = VariableFeatures(object = combined_seurat))
combined_seurat <- RunUMAP(combined_seurat, dims = 1:30)
combined_seurat <- RunTSNE(combined_seurat, dims = 1:30)

# Clustering
cat("Performing clustering...\n")
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:30)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

# Save cleaned and integrated data
cat("\n=== SAVING CLEANED DATA ===\n")

# Save Seurat object
saveRDS(combined_seurat, "results/cleaned_data/GSE222520_cleaned_seurat.rds")

# Save cell metadata
cell_metadata <- combined_seurat@meta.data
write.csv(cell_metadata, "data/processed/cell_metadata/GSE222520_cell_metadata.csv", row.names = TRUE)

# Save gene metadata
gene_metadata <- data.frame(
  gene_id = rownames(combined_seurat),
  gene_name = rownames(combined_seurat),
  stringsAsFactors = FALSE
)
write.csv(gene_metadata, "data/processed/gene_metadata/GSE222520_gene_metadata.csv", row.names = FALSE)

# Save expression matrix (sparse format for efficiency)
expression_matrix <- GetAssayData(combined_seurat, slot = "counts")
saveRDS(expression_matrix, "data/processed/expression_matrices/GSE222520_expression_matrix.rds")

# Create summary plots
cat("Creating summary plots...\n")

# UMAP colored by sample
p1 <- DimPlot(combined_seurat, reduction = "umap", group.by = "sample_id", pt.size = 0.1) +
  ggtitle("UMAP by Sample")

# UMAP colored by group
p2 <- DimPlot(combined_seurat, reduction = "umap", group.by = "group", pt.size = 0.1) +
  ggtitle("UMAP by Group")

# UMAP colored by diagnosis
p3 <- DimPlot(combined_seurat, reduction = "umap", group.by = "diagnosis", pt.size = 0.1) +
  ggtitle("UMAP by Diagnosis")

# UMAP colored by genotype
p4 <- DimPlot(combined_seurat, reduction = "umap", group.by = "genotype", pt.size = 0.1) +
  ggtitle("UMAP by Genotype")

# Save plots
ggsave("results/qc_plots/umap_by_sample.pdf", p1, width = 10, height = 8)
ggsave("results/qc_plots/umap_by_group.pdf", p2, width = 10, height = 8)
ggsave("results/qc_plots/umap_by_diagnosis.pdf", p3, width = 10, height = 8)
ggsave("results/qc_plots/umap_by_genotype.pdf", p4, width = 10, height = 8)

# Create comprehensive summary
cat("\n=== FINAL SUMMARY ===\n")
cat("Total samples processed:", length(seurat_objects), "\n")
cat("Total cells after filtering:", ncol(combined_seurat), "\n")
cat("Total genes:", nrow(combined_seurat), "\n")
cat("Variable features identified:", length(VariableFeatures(combined_seurat)), "\n")
cat("Clusters identified:", length(unique(combined_seurat$seurat_clusters)), "\n")

# Sample distribution
cat("\nSample distribution:\n")
sample_counts <- table(combined_seurat$sample_id)
for (sample in names(sample_counts)) {
  cat("  ", sample, ":", sample_counts[sample], "cells\n")
}

# Group distribution
cat("\nGroup distribution:\n")
group_counts <- table(combined_seurat$group)
for (group in names(group_counts)) {
  cat("  ", group, ":", group_counts[group], "cells\n")
}

# Save processing summary
summary_data <- list(
  total_samples = length(seurat_objects),
  total_cells = ncol(combined_seurat),
  total_genes = nrow(combined_seurat),
  variable_features = length(VariableFeatures(combined_seurat)),
  clusters = length(unique(combined_seurat$seurat_clusters)),
  sample_distribution = sample_counts,
  group_distribution = group_counts,
  processing_date = Sys.Date()
)

saveRDS(summary_data, "results/cleaned_data/processing_summary.rds")

cat("\n=== DATA CLEANING COMPLETE ===\n")
cat("Cleaned data saved to: results/cleaned_data/\n")
cat("QC plots saved to: results/qc_plots/\n")
cat("Processed data saved to: data/processed/\n")
