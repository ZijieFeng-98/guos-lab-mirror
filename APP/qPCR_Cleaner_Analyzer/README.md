# Immunology Analysis Project

## Overview
This project performs comprehensive correlation analysis between pathway genes and immune markers using TCGA cancer genomics data across multiple cancer types.

## 🚀 Quick Start

### 1. Install Dependencies
```bash
Rscript scripts/utilities/install_packages.R
```

### 2. Run Analysis
```bash
# Run GBM analysis
Rscript scripts/analysis/TCGA_GBM_Correlation_Analysis_Expanded_Immune.R

# Run other cancer types
Rscript scripts/analysis/TCGA_BRCA_Correlation_Analysis.R
Rscript scripts/analysis/TCGA_OV_Correlation_Analysis.R
Rscript scripts/analysis/TCGA_PAAD_Correlation_Analysis.R
```

### 3. View Results
- **Figures**: `results/figures/` - Heatmaps and visualizations
- **Tables**: `results/tables/` - Correlation matrices

## 📁 Project Structure

```
Immunology_analysis/
├── scripts/
│   ├── analysis/          # Main analysis scripts
│   ├── data_download/     # Data download utilities
│   └── utilities/         # Setup and utility scripts
├── data/
│   ├── raw/              # Raw TCGA data files
│   └── processed/        # Processed data
├── results/
│   ├── figures/          # Generated plots
│   └── tables/           # Statistical results
└── docs/                 # Documentation
```

📖 **Detailed structure**: See [docs/PROJECT_STRUCTURE.md](docs/PROJECT_STRUCTURE.md)

## 🧬 Cancer Types Analyzed

| Cancer Type | Code | Description |
|-------------|------|-------------|
| GBM | Glioblastoma Multiforme | Brain cancer |
| BRCA | Breast Cancer | Breast cancer |
| OV | Ovarian Cancer | Ovarian cancer |
| PAAD | Pancreatic Cancer | Pancreatic cancer |

## 🔬 Analysis Components

### Pathway Genes (19 genes)
- **Cholesterol/Lipid Metabolism**: FLVCR1, FLVCR2, TXNDC16, SOAT1, SOAT2, SCAP, SREBF1, LRP1, LRP8, LDLR, APOA1, APOB, APOC1, APOD, APOE, DGAT1, DGAT2
- **Innate Immunity**: TMEM173, MB21D1

### Immune Markers (50+ genes)
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

## 📊 Analysis Output

### Generated Files
- **Correlation Heatmaps**: PNG files showing gene-gene correlations
- **Correlation Matrices**: CSV files with numerical correlation values
- **Statistical Results**: Comprehensive correlation analysis

### File Naming Convention
- `Spearman_Correlation_Heatmap_{CANCER_TYPE}_All_Immune.png`
- `Spearman_Correlation_Matrix_{CANCER_TYPE}_All_Immune.csv`

## 🛠️ Requirements

- **R version**: 3.6 or higher
- **Key packages**: TCGAbiolinks, SummarizedExperiment, pheatmap
- **Disk space**: ~2-5 GB for TCGA data
- **Memory**: 8GB+ RAM recommended

## 📚 Documentation

- **[Project Structure](docs/PROJECT_STRUCTURE.md)**: Detailed directory organization
- **[Cleanup Summary](docs/CLEANUP_SUMMARY.md)**: Project cleanup and space optimization
- **[Original README](docs/README.md)**: Original project documentation
- **[Manual Download Guide](docs/manual_download_guide.md)**: Offline data download instructions

## 🔧 Troubleshooting

### GDC Server Issues
If you encounter "GDC server down" errors:
1. Use offline analysis scripts: `*_Offline.R`
2. Check server status: `Rscript scripts/utilities/check_gdc_status.R`
3. Follow manual download guide in `docs/`

### Memory Issues
- Ensure sufficient RAM (8GB+ recommended)
- Close other applications during analysis
- Consider running on a high-memory machine

## 📈 Results Interpretation

- **Positive correlations** (red in heatmap): Genes that tend to be expressed together
- **Negative correlations** (blue in heatmap): Genes with opposite expression patterns
- **No correlation** (white in heatmap): Genes with independent expression

## 🤝 Contributing

1. Follow the organized project structure
2. Place new scripts in appropriate directories
3. Update documentation as needed
4. Use consistent naming conventions

## 📄 License

This project is for research purposes. Please cite original sources when using results.
