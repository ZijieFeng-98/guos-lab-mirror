# GSE222520 Data Summary Report

## Dataset Overview
**GEO Accession**: GSE222520
**Download Date**: 2025-08-21 13:15:28
**Total Samples**: 18
**Total Cells**: 142,737
**Total Genes**: 33,538
**Total UMIs**: 832,218,956

## Sample Distribution
| Sample Type | Count | Description |
|-------------|-------|-------------|
| NGB | 3 | Normal brain tissue |
| IMP | 3 | Primary glioma |
| IMR | 6 | Recurrent glioma |
| IWP | 4 | Primary glioma with treatment |
| IWR | 3 | Recurrent glioma with treatment |

## Quality Metrics
| Metric | Value |
|--------|-------|
| Mean UMIs per cell | 5830.4 |
| Median UMIs per cell | 4537 |
| Mean genes per cell | 7929.8 |

## Data Files
### Expression Data
- `GSE222520_combined_expression_matrix.rds`: Sparse matrix format
- `GSE222520_combined_expression_matrix_test.csv`: CSV format (first 1000 genes)

### Metadata
- `GSE222520_cell_metadata.csv`: Cell-level information
- `GSE222520_gene_metadata.csv`: Gene-level information
- `sample_info.csv`: Sample-level information

## Processing Notes
- Data downloaded from GEO as processed 10x Genomics format
- Successfully loaded 18 out of 21 samples (3 had corrupted files)
- All samples have consistent gene set (33,538 genes)
- Cell barcodes prefixed with sample names to avoid conflicts

## Next Steps
1. Quality control and filtering
2. Normalization and scaling
3. Feature selection
4. Dimensionality reduction
5. Clustering analysis
6. Cell type annotation
7. Differential expression analysis
