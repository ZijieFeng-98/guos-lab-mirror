# GSE222520 Single-Cell RNA-seq Analysis Project

## Project Overview
This project contains processed single-cell RNA-seq data from GSE222520, which analyzes the glioma immune microenvironment.

## Dataset Information
- **GEO Accession**: GSE222520
- **Title**: Single-cell RNA-seq analysis of glioma immune microenvironment
- **Data Type**: Processed scRNA-seq (10x Genomics)
- **Total Cells**: 142,737
- **Total Genes**: 33,538
- **Total UMIs**: 832,218,956

## Sample Types
- **NGB**: Normal brain tissue (3 samples)
- **IMP**: Primary glioma (3 samples)
- **IMR**: Recurrent glioma (6 samples)
- **IWP**: Primary glioma with treatment (4 samples)
- **IWR**: Recurrent glioma with treatment (3 samples)

## Project Structure
```
Sc_Immune_Analysis/
├── data/
│   ├── raw/                    # Raw downloaded files
│   ├── processed/              # Processed expression matrices
│   ├── metadata/               # Cell and sample metadata
│   └── annotations/            # Gene annotations
├── scripts/
│   ├── download/               # Data download scripts
│   └── analysis/               # Analysis scripts
├── results/
│   ├── qc/                     # Quality control results
│   └── clustering/             # Clustering results
├── figures/
│   ├── qc/                     # QC plots
│   └── clustering/             # Clustering plots
└── docs/
    └── data_description/       # Data documentation
```

## Key Files
- `data/processed/GSE222520_combined_expression_matrix.rds`: Combined expression matrix
- `data/processed/GSE222520_cell_metadata.csv`: Cell-level metadata
- `data/processed/GSE222520_gene_metadata.csv`: Gene-level metadata
- `data/processed/sample_info.csv`: Sample-level information

## Data Processing Steps
1. Downloaded processed data from GEO (GSE222520)
2. Extracted tar.gz files containing 10x Genomics data
3. Loaded and combined expression matrices from 18 samples
4. Created comprehensive metadata files

## Notes
- Data is already processed and batch-corrected (Harmony applied)
- Cell type annotations are available
- Reference signatures (GlioTIME-36) are included
- Ready for secondary analysis (differential expression, pathway analysis, etc.)

## Scripts Used
- `download_processed_gse222520.R`: Initial download and organization
- `extract_and_organize_data.R`: Extract tar.gz files
- `load_and_combine_matrices_fixed.R`: Load and combine expression matrices

Generated on:
2025-08-21 13:15:28
