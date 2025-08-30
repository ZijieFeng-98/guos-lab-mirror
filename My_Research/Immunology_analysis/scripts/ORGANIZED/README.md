# IMMUNOLOGY ANALYSIS SCRIPTS - ORGANIZED BY PURPOSE

## Overview
This directory contains all analysis scripts organized by their specific purpose and function.

## Directory Structure

```
scripts/ORGANIZED/
├── ANALYSIS/              # Main correlation analysis scripts
│   ├── TCGA_GBM_Correlation_Analysis_Offline.R
│   ├── TCGA_PAAD_Correlation_Analysis_Offline.R
│   ├── TCGA_OV_Correlation_Analysis_Offline.R
│   ├── TCGA_BRCA_Correlation_Analysis_Offline.R
│   ├── TCGA_GBM_Correlation_Analysis.R
│   ├── TCGA_PAAD_Correlation_Analysis.R
│   ├── TCGA_OV_Correlation_Analysis.R
│   ├── TCGA_BRCA_Correlation_Analysis.R
│   └── TCGA_GBM_Correlation_Analysis_Expanded_Immune.R
├── DATA_DOWNLOAD/         # Data download scripts
│   ├── download_tcga_data.R
│   ├── download_paad_data.R
│   ├── download_ov_data.R
│   ├── download_ov_simple.R
│   └── download_clinical_data_simple.R
└── UTILITIES/             # Utility and check scripts
    ├── check_gdc_status.R
    ├── check_ov_genes.R
    ├── check_ov_gene_names.R
    └── install_packages.R
```

## Script Categories

### ANALYSIS Scripts (8 files)
**Purpose**: Perform correlation analysis between target genes and immune markers
- **Offline versions**: Use pre-downloaded data (.rds files)
- **Online versions**: Download data from TCGA directly
- **Expanded Immune**: Additional immune markers included

### DATA_DOWNLOAD Scripts (5 files)
**Purpose**: Download and prepare TCGA data
- Download RNA-seq data for different cancer types
- Download clinical data
- Simple and comprehensive versions available

### UTILITIES Scripts (4 files)
**Purpose**: Support and maintenance functions
- Check GDC status and connectivity
- Verify gene names and data integrity
- Install required packages

## Usage

### For Analysis:
```bash
# Run GBM analysis with offline data
Rscript scripts/ORGANIZED/ANALYSIS/TCGA_GBM_Correlation_Analysis_Offline.R

# Run PAAD analysis with online data
Rscript scripts/ORGANIZED/ANALYSIS/TCGA_PAAD_Correlation_Analysis.R
```

### For Data Download:
```bash
# Download GBM data
Rscript scripts/ORGANIZED/DATA_DOWNLOAD/download_tcga_data.R
```

### For Utilities:
```bash
# Check GDC status
Rscript scripts/ORGANIZED/UTILITIES/check_gdc_status.R

# Install packages
Rscript scripts/ORGANIZED/UTILITIES/install_packages.R
```

## File Organization
All scripts are organized by their specific purpose for easy navigation and maintenance.
