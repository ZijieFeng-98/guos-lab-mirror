# TCGA GBM Correlation Analysis: Expanded Immune Marker List

## Overview
This project performs correlation analysis between pathway genes and an expanded immune marker list using TCGA-GBM (Glioblastoma) transcriptome data. The analysis is based on the merged immune pathway markers document for glioma research.

## Files
- `TCGA_GBM_Correlation_Analysis_Expanded_Immune.R`: Main R script for the correlation analysis
- `README.md`: This documentation file

## Analysis Components

### Pathway Genes (19 genes)
- **Cholesterol/Lipid Metabolism**: FLVCR1, FLVCR2, TXNDC16, SOAT1, SOAT2, SCAP, SREBF1, LRP1, LRP8, LDLR, APOA1, APOB, APOC1, APOD, APOE, DGAT1, DGAT2
- **Innate Immunity**: TMEM173, MB21D1

### Immune Markers (50+ genes)
Organized into functional categories:
- **Cytokines/Cytolytic**: IFNA1, IFNB1, IL1B, IL6, TNF, GZMA, GZMB, PRF1, CXCL9, CXCL10, CCL2, CCL4, CCL5
- **IFN-related**: STAT1, IRF7, ISG15, IFIT1, IFI6, MX1, OAS1
- **Immune Checkpoints**: CD274, HAVCR2, IDO1, LGALS9
- **MHC/Antigen Presentation**: HLA-DRA, HLA-DRB1, CD74, CD83, CTSB, CTSS
- **T Cell Markers**: CD3D, CD4, CD8A
- **B Cell Markers**: CD19, CD79A
- **Macrophage/Monocyte**: CD14, CD68, CD163, MRC1, ITGAM, CCR2, FCN1, LYZ
- **Microglia**: P2RY12, TMEM119, CX3CR1, TREM2, GPR34
- **Dendritic Cell**: CD1C, ITGAX, CLEC9A, CD86
- **NK Cell**: NCAM1, NKG7, KLRD1, GNLY

## Prerequisites

### Required R Packages
```r
install.packages(c("TCGAbiolinks", "SummarizedExperiment", "pheatmap"))
```

### System Requirements
- R version 3.6 or higher
- Sufficient disk space for TCGA data download (~2-5 GB)
- Internet connection for data download

## Usage

1. **Set up R environment**:
   ```r
   # Install required packages if not already installed
   if (!require(TCGAbiolinks)) install.packages("TCGAbiolinks")
   if (!require(SummarizedExperiment)) install.packages("SummarizedExperiment")
   if (!require(pheatmap)) install.packages("pheatmap")
   ```

2. **Run the analysis**:
   ```r
   source("TCGA_GBM_Correlation_Analysis_Expanded_Immune.R")
   ```

## Output Files

The script generates the following outputs in the `Processed_Data/` directory:

1. **Spearman_Correlation_Heatmap_All_Immune.png**: High-resolution heatmap visualization
2. **Spearman_Correlation_Matrix_All_Immune.csv**: Raw correlation matrix for further analysis

## Analysis Steps

1. **Data Download**: Queries and downloads TCGA-GBM transcriptome data
2. **Sample Filtering**: Keeps only primary tumor samples (TP)
3. **Gene Filtering**: Extracts expression data for pathway and immune marker genes
4. **Correlation Analysis**: Computes Spearman correlations between all gene pairs
5. **Visualization**: Creates clustered heatmap
6. **Data Export**: Saves results for further analysis

## Notes

- The analysis focuses on primary tumor samples only
- Gene names are cleaned to remove duplicates and NA values
- The heatmap uses hierarchical clustering for both rows and columns
- Correlation values range from -1 to +1 (negative to positive correlation)

## Troubleshooting

- **Memory issues**: Consider running on a machine with sufficient RAM (8GB+ recommended)
- **Download errors**: Check internet connection and TCGA server status
- **Package installation**: Ensure R is up to date and packages are compatible

### GDC Server Issues

If you encounter "GDC server down" errors:

1. **Check server status**: Visit https://portal.gdc.cancer.gov/
2. **Use offline mode**: Run `source('TCGA_GBM_Correlation_Analysis_Offline.R')` if you have previously downloaded data
3. **Download when available**: Use `source('download_tcga_data.R')` when the server is back online
4. **Manual download**: Follow instructions in the offline script for manual data download

### Alternative Data Sources

If GDC remains unavailable, consider:
- Using previously downloaded TCGA data
- Alternative cancer genomics databases (ICGC, GEO)
- Local institutional data repositories

## Citation
This analysis is based on the merged immune pathway markers document for glioma research. Please cite the original sources when using these results.
