# Data Organization

This directory contains all data files for the Immunology Analysis project.

## Directory Structure

### `raw/`
Contains all raw data files:

- **`*.rds`**: R data files containing processed TCGA data
  - `brca_data.rds`: Breast cancer data
  - `gbm_data.rds`: Glioblastoma data (original)
  - `gbm_data_new.rds`: Glioblastoma data (updated)
  - `ov_data.rds`: Ovarian cancer data
  - `paad_data.rds`: Pancreatic adenocarcinoma data

- **`rna_seq_counts/`**: RNA-seq gene count files organized by sample UUIDs
  - Each UUID directory contains a `.rna_seq.augmented_star_gene_counts.tsv` file

- **`archives/`**: Compressed data archives
  - `Mon_Aug_18_10_36_31_2025_0.tar.gz`
  - `Mon_Aug_18_10_36_31_2025_1.tar.gz`

### `processed/`
Contains processed/derived data files (currently empty)

## Data Sources

- **TCGA Data**: Downloaded from GDC (Genomic Data Commons)
- **RNA-seq Counts**: Gene expression quantification data
- **Cancer Types**: BRCA (Breast), GBM (Glioblastoma), OV (Ovarian), PAAD (Pancreatic)

## Notes

- The `GDCdata/` directory in the project root contains the original TCGA downloads
- RNA-seq count files were previously scattered in UUID-named directories in the root
- Duplicate .rds files have been consolidated in `data/raw/`
