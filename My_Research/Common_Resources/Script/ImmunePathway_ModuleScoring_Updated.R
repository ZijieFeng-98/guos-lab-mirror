
# -------------------------------
# scRNA-seq Immune Pathway Scoring (Based on Wang et al. + GSE222520 + Custom Immune Targets)
# Author: Guo Lab | Adapted for Zijie Feng
# -------------------------------

# 📦 Load Libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# 📂 Load Processed Seurat Object
# Replace this with the correct path to your GSE163120 dataset
# Example:
# seurat_obj <- readRDS("path_to_seurat_object.rds")

# ------------------------------------
# 🧬 Define Gene Sets for Module Scoring
# ------------------------------------

# APOE Pathway
apoe_genes <- c("APOE", "APOC1", "LRP1")

# SREBP / Lipid Metabolism
srebp_genes <- c("SCD", "FASN", "HMGCR", "LDLR", "INSIG1", "ACACA", "ACLY")

# Interferon Response
ifn_genes <- c("STAT1", "CXCL10", "IFIT1", "ISG15", "IRF7", "IFI6")

# Hypoxia + HILPDA
hypoxia_genes <- c("HILPDA", "VEGFA", "CA9", "LDHA", "ENO1")

# Hexosamine Pathway / GFPT1
hex_genes <- c("GFPT1", "OGT", "MGEA5", "HSP90AA1", "NFE2L2")

# Cytokines / Chemokines
cytokine_genes <- c("IL1B", "CCL2", "CXCL9", "CXCL10", "TNF")

# Phagocytosis / Antigen Presentation
phagocytosis_genes <- c("CD68", "CD74", "CTSB", "CTSS", "HLA-DRA", "HLA-DRB1")

# Immunosuppression / Checkpoint
checkpoint_genes <- c("CD274", "LGALS9", "HAVCR2", "IDO1")

# ------------------------------------
# Additional Modules: GSE222520 Pathway Targets
# ------------------------------------

# Lipid Transport & Uptake
lipid_uptake_genes <- c("FLVCR1", "FLVCR2", "TXNDC16", "LRP1", "LRP8", "LDLR")

# Cholesterol Biosynthesis (SREBP Pathway)
srebp_pathway_genes <- c("SCAP", "SREBF1", "SOAT1", "SOAT2", "DGAT1", "DGAT2")

# Apolipoproteins
apolipoprotein_genes <- c("APOA1", "APOB", "APOC1", "APOD", "APOE")

# STING Pathway (Innate Immune Signaling)
sting_genes <- c("TMEM173")

# ------------------------------------
# 🧮 Add Module Scores
# ------------------------------------
seurat_obj <- AddModuleScore(seurat_obj, features = list(apoe_genes), name = "APOE_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(srebp_genes), name = "SREBP_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(ifn_genes), name = "IFN_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(hypoxia_genes), name = "HILPDA_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(hex_genes), name = "HEX_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(cytokine_genes), name = "CYTOKINE_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(phagocytosis_genes), name = "PHAGO_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(checkpoint_genes), name = "CHECK_Score")

seurat_obj <- AddModuleScore(seurat_obj, features = list(lipid_uptake_genes), name = "LipidUptake_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(srebp_pathway_genes), name = "SREBP_Pathway_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(apolipoprotein_genes), name = "APO_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(sting_genes), name = "STING_Score")

# ------------------------------------
# 📊 Visualization Examples
# ------------------------------------
VlnPlot(seurat_obj, features = c("APOE_Score1", "SREBP_Score1", "IFN_Score1"), group.by = "cell_type")
FeaturePlot(seurat_obj, features = c("HILPDA_Score1", "HEX_Score1", "LipidUptake_Score1"), split.by = "condition")
VlnPlot(seurat_obj, features = c("SREBP_Pathway_Score1", "STING_Score1"), group.by = "condition")

# Save Seurat object (optional)
# saveRDS(seurat_obj, "annotated_seurat_with_scores.rds")
