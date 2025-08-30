# Immunology Analysis Project Structure

## Overview
This project performs correlation analysis between pathway genes and immune markers using TCGA cancer genomics data across multiple cancer types (GBM, BRCA, OV, PAAD).

## Directory Structure

```
Immunology_analysis/
├── scripts/
│   ├── analysis/                    # Main analysis scripts
│   │   ├── TCGA_BRCA_Correlation_Analysis.R
│   │   ├── TCGA_BRCA_Correlation_Analysis_Offline.R
│   │   ├── TCGA_GBM_Correlation_Analysis_Expanded_Immune.R
│   │   ├── TCGA_GBM_Correlation_Analysis_Offline.R
│   │   ├── TCGA_OV_Correlation_Analysis.R
│   │   ├── TCGA_OV_Correlation_Analysis_Offline.R
│   │   ├── TCGA_PAAD_Correlation_Analysis.R
│   │   └── TCGA_PAAD_Correlation_Analysis_Offline.R
│   ├── data_download/               # Data download scripts
│   │   ├── download_tcga_data.R
│   │   ├── download_ov_data.R
│   │   ├── download_ov_simple.R
│   │   └── download_paad_data.R
│   └── utilities/                   # Utility and setup scripts
│       ├── install_packages.R
│       ├── check_gdc_status.R
│       ├── check_ov_gene_names.R
│       └── check_ov_genes.R
├── data/
│   ├── raw/                         # Raw data files
│   │   ├── brca_data.rds
│   │   ├── gbm_data.rds
│   │   ├── ov_data.rds
│   │   └── paad_data.rds
│   └── processed/                   # Processed data (if any)
├── results/
│   ├── figures/                     # Generated plots and visualizations
│   │   ├── Spearman_Correlation_Heatmap_All_Immune.png
│   │   ├── Spearman_Correlation_Heatmap_BRCA_All_Immune.png
│   │   ├── Spearman_Correlation_Heatmap_OV_All_Immune.png
│   │   └── Spearman_Correlation_Heatmap_PAAD_All_Immune_Offline.png
│   └── tables/                      # Generated data tables
│       ├── Spearman_Correlation_Matrix_All_Immune.csv
│       ├── Spearman_Correlation_Matrix_BRCA_All_Immune.csv
│       ├── Spearman_Correlation_Matrix_OV_All_Immune.csv
│       └── Spearman_Correlation_Matrix_PAAD_All_Immune_Offline.csv
├── docs/                            # Documentation
│   ├── README.md
│   ├── PROJECT_STRUCTURE.md
│   └── manual_download_guide.md
├── GDCdata/                         # TCGA data directory
├── Immunology_analysis.Rproj        # RStudio project file
├── MANIFEST.txt                     # Project manifest
├── run_analysis.sh                  # Shell script for running analysis
└── .RData                           # R workspace file
```

## File Descriptions

### Analysis Scripts (`scripts/analysis/`)
- **TCGA_*_Correlation_Analysis.R**: Main analysis scripts for each cancer type
- **TCGA_*_Correlation_Analysis_Offline.R**: Offline versions for when GDC server is down

### Data Download Scripts (`scripts/data_download/`)
- **download_tcga_data.R**: General TCGA data download utility
- **download_*_data.R**: Cancer-specific data download scripts

### Utility Scripts (`scripts/utilities/`)
- **install_packages.R**: Install and check required R packages
- **check_*.R**: Various utility scripts for data validation and status checking

### Data Files (`data/`)
- **raw/**: Contains the large `.rds` files with TCGA transcriptome data
- **processed/**: For any intermediate processed data files

### Results (`results/`)
- **figures/**: Generated heatmaps and visualizations
- **tables/**: Correlation matrices and statistical results

## Usage

### Running Analysis
```bash
# Run GBM analysis
Rscript scripts/analysis/TCGA_GBM_Correlation_Analysis_Expanded_Immune.R

# Run other cancer types
Rscript scripts/analysis/TCGA_BRCA_Correlation_Analysis.R
Rscript scripts/analysis/TCGA_OV_Correlation_Analysis.R
Rscript scripts/analysis/TCGA_PAAD_Correlation_Analysis.R
```

### Installing Dependencies
```bash
Rscript scripts/utilities/install_packages.R
```

### Downloading Data
```bash
Rscript scripts/data_download/download_tcga_data.R
```

## Cancer Types Analyzed
1. **GBM** (Glioblastoma Multiforme)
2. **BRCA** (Breast Cancer)
3. **OV** (Ovarian Cancer)
4. **PAAD** (Pancreatic Cancer)

## Analysis Components
- **Pathway Genes**: 19 genes related to cholesterol/lipid metabolism and innate immunity
- **Immune Markers**: 50+ genes across multiple immune cell types and functions
- **Correlation Analysis**: Spearman correlations between all gene pairs
- **Visualization**: Hierarchical clustering heatmaps

