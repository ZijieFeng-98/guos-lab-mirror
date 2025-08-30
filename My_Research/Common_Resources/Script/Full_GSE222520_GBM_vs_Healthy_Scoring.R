
# -------------------------------
# Full scRNA-seq Workflow for GBM vs Healthy Immune Analysis
# Based on GSE222520 | For use in Seurat (R)
# Author: Guo Lab | Adapted for Zijie Feng
# -------------------------------

# 📦 Load Libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# -------------------------------
# Step 1: Load Processed Data
# -------------------------------
# Replace with your actual directory paths
sc_data <- Read10X(data.dir = "path_to_filtered_matrix/")
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "GSE222520", min.cells = 3, min.features = 200)

# Load and merge metadata if available
# metadata <- read.csv("GSE222520_metadata.csv")
# seurat_obj <- AddMetaData(seurat_obj, metadata)

# -------------------------------
# Step 2: Preprocessing & QC
# -------------------------------
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# -------------------------------
# Step 3: (Optional) Integration if analyzing per-sample
# -------------------------------
# sample.list <- SplitObject(seurat_obj, split.by = "sample")
# sample.list <- lapply(sample.list, SCTransform)
# features <- SelectIntegrationFeatures(sample.list, nfeatures = 3000)
# anchors <- FindIntegrationAnchors(sample.list, anchor.features = features)
# seurat_obj <- IntegrateData(anchorset = anchors)

# -------------------------------
# Step 4: PCA, Clustering, UMAP
# -------------------------------
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Visualize by condition or diagnosis
DimPlot(seurat_obj, group.by = "diagnosis")  # Adjust column name as needed

# -------------------------------
# Step 5: Immune Cell Annotation
# -------------------------------
# Visualize known markers
FeaturePlot(seurat_obj, features = c("CD14", "P2RY12", "CD3D", "CD8A"))
# After visual review, annotate clusters manually if needed

# -------------------------------
# Step 6: Cell Type Proportion Comparison
# -------------------------------
cell_props <- seurat_obj@meta.data %>%
  group_by(diagnosis, seurat_clusters) %>%
  summarise(n = n()) %>%
  group_by(diagnosis) %>%
  mutate(freq = n / sum(n))

ggplot(cell_props, aes(x = diagnosis, y = freq, fill = factor(seurat_clusters))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(fill = "Cluster") +
  theme_minimal()

# -------------------------------
# Step 7: Pathway Scoring (FLVCR1–DGAT2 Range + Others)
# -------------------------------

# Grouped gene modules
srebp_genes <- c("SCAP", "SREBF1", "SOAT1", "SOAT2", "DGAT1", "DGAT2")
apolipo_genes <- c("APOA1", "APOB", "APOC1", "APOD", "APOE", "LDLR", "LRP1", "LRP8")
sting_cgas_genes <- c("TMEM173", "MB21D1")
heme_redox_genes <- c("FLVCR1", "FLVCR2", "TXNDC16")

# Add module scores
seurat_obj <- AddModuleScore(seurat_obj, features = list(srebp_genes), name = "SREBP_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(apolipo_genes), name = "APO_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(sting_cgas_genes), name = "STING_cGAS_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(heme_redox_genes), name = "HEME_Score")

# -------------------------------
# Step 8: Visualization
# -------------------------------
# Violin plots to compare pathway activity between GBM vs Healthy
VlnPlot(seurat_obj, features = c("SREBP_Score1", "APO_Score1", "STING_cGAS_Score1", "HEME_Score1"), group.by = "diagnosis")

# Optional UMAPs
FeaturePlot(seurat_obj, features = c("APOE", "SREBF1", "TMEM173", "FLVCR1"), split.by = "diagnosis")
