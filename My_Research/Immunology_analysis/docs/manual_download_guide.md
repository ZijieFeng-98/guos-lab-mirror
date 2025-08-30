# Manual TCGA Data Download Guide

## When GDC Server is Down

This guide provides step-by-step instructions for manually downloading TCGA-GBM data when the automated tools are unavailable.

## Method 1: GDC Data Portal (Web Interface)

### Step 1: Access the Portal
1. Visit: https://portal.gdc.cancer.gov/
2. Click on "Repository" in the top navigation
3. Click on "Files"

### Step 2: Apply Filters
Set the following filters in the left sidebar:

**Project:**
- Select "TCGA-GBM"

**Data Category:**
- Select "Transcriptome Profiling"

**Data Type:**
- Select "Gene Expression Quantification"

**Workflow Type:**
- Select "STAR - Counts"

**Sample Type:**
- Select "Primary Tumor"

### Step 3: Download Files
1. Click "Add All Files to Cart" (or select specific files)
2. Click the cart icon in the top right
3. Click "Download Manifest"
4. Save the manifest file (e.g., `gdc_manifest.txt`)

### Step 4: Use GDC Data Transfer Tool
1. Download GDC Data Transfer Tool from: https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
2. Install the tool
3. Open terminal/command prompt
4. Navigate to your project directory
5. Run: `gdc-client download -m gdc_manifest.txt`

## Method 2: Using R with Manual Download

### Step 1: Download Individual Files
If the portal allows individual file downloads:

1. Download the count files (`.txt` or `.tsv` format)
2. Save them in a folder named `TCGA_GBM_Data`

### Step 2: Create R Script to Process Files
```r
# Manual data processing script
library(SummarizedExperiment)

# Set path to downloaded files
data_path <- "TCGA_GBM_Data/"

# List all count files
count_files <- list.files(data_path, pattern = "\\.txt$|\\.tsv$", full.names = TRUE)

# Read and combine data
expression_data <- list()
sample_ids <- c()

for (file in count_files) {
  # Extract sample ID from filename
  sample_id <- gsub(".*/([^/]+)\\..*", "\\1", file)
  sample_ids <- c(sample_ids, sample_id)
  
  # Read count data
  counts <- read.table(file, header = TRUE, sep = "\t")
  expression_data[[sample_id]] <- counts
}

# Combine all samples
# (This is a simplified version - actual implementation depends on file format)
```

## Method 3: Alternative Data Sources

### Option A: GEO Database
1. Visit: https://www.ncbi.nlm.nih.gov/geo/
2. Search for "TCGA GBM" or "glioblastoma"
3. Look for datasets with gene expression data
4. Download processed data files

### Option B: ICGC Data Portal
1. Visit: https://dcc.icgc.org/
2. Search for glioblastoma/glioma datasets
3. Download available gene expression data

### Option C: Local Institutional Data
- Check if your institution has local copies of TCGA data
- Contact your bioinformatics core facility
- Check shared research drives

## Method 4: Using Sample Data for Testing

If you need to test your analysis pipeline immediately:

```r
# Run the sample data generator
source('create_sample_data.R')

# This will create synthetic data that mimics TCGA structure
# You can then test your analysis pipeline
source('TCGA_GBM_Correlation_Analysis_Offline.R')
```

## Converting Downloaded Data to RDS Format

Once you have the raw data files, convert them to the expected format:

```r
# Example conversion script
library(SummarizedExperiment)

# Assuming you have expression matrix and metadata
expression_matrix <- your_expression_data
sample_metadata <- your_sample_metadata
gene_metadata <- your_gene_metadata

# Create SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(counts = expression_matrix),
  rowData = gene_metadata,
  colData = sample_metadata
)

# Save as RDS file
saveRDS(se, file = "gbm_data.rds")
```

## Troubleshooting

### Common Issues:
1. **Large file downloads**: Use download managers or command-line tools
2. **File format issues**: Check file headers and separators
3. **Memory limitations**: Process files in batches
4. **Network timeouts**: Use resume-capable download tools

### File Format Requirements:
- Gene expression data should be in matrix format (genes × samples)
- Gene names should be in standard format (e.g., gene symbols)
- Sample IDs should be consistent
- Missing values should be handled appropriately

## Next Steps

After obtaining the data:
1. Save as `gbm_data.rds` in your project directory
2. Run the offline analysis: `source('TCGA_GBM_Correlation_Analysis_Offline.R')`
3. Verify results and adjust parameters as needed

## Notes

- Manual downloads may take significant time depending on file sizes
- Always verify data integrity after download
- Keep backup copies of downloaded files
- Consider using sample data for initial testing and development
