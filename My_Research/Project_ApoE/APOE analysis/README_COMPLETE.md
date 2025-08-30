# APOE Isoform Analysis - Complete Workflow

## Overview
This repository contains a comprehensive APOE (Apolipoprotein E) isoform analysis workflow that can analyze genomic data to determine APOE ε2, ε3, and ε4 allele frequencies. The workflow is designed to work with 1000 Genomes Project data and can be adapted for TCGA cancer genomics studies.

## Quick Start

### Prerequisites
- R (version 4.0 or higher)
- Required R packages (see Installation section)

### Installation
1. Install required R packages:
```r
install.packages(c("vcfR", "dplyr", "tidyr", "ggplot2", "readr", "stringr", "purrr", "forcats"))
```

2. Clone or download this repository

### Running the Analysis
```bash
Rscript run_apoe_analysis.R
```

This will:
- Download 1000 Genomes data (if not already present)
- Extract APOE SNPs (rs429358, rs7412)
- Assign APOE isoforms based on haplotype logic
- Generate comprehensive results and visualizations

## Files Overview

### Core Analysis Scripts
- `apoe_analysis_complete.R` - Complete analysis workflow with full documentation
- `run_apoe_analysis.R` - Simple runner script
- `run_apoe_simple.R` - Simplified version using pre-extracted data

### Documentation
- `APOE_ANALYSIS_SUMMARY.md` - Comprehensive summary of logic, results, and methodology
- `README_COMPLETE.md` - This file

### Output Files (Generated)
- `output/apoe_genotypes.csv` - Individual-level genotype data
- `output/haplotype_frequencies.csv` - Overall haplotype frequencies
- `output/genotype_frequencies.csv` - Overall genotype frequencies
- `output/superpopulation_summary.csv` - Population-stratified results
- `output/*.png` - Visualization plots

## APOE Logic Summary

### Key SNPs
- **rs429358**: chr19:45411941 (T→C)
- **rs7412**: chr19:45412079 (C→T)

### Haplotype Assignment
```
rs429358 | rs7412 | APOE Isoform
---------|--------|-------------
T        | T      | ε2
T        | C      | ε3
C        | C      | ε4
C        | T      | ε1 (rare)
```

### Genotype Assignment
- Individual genotypes combine two haplotypes (e.g., ε3/ε4)
- Sorted alphabetically (ε2/ε4 not ε4/ε2)

## Key Results from 1000 Genomes

### Overall Frequencies
- **ε3**: 77.5% (most common)
- **ε4**: 15.0%
- **ε2**: 7.5%
- **ε1**: 0.02% (rare)

### Population Differences
- **African populations**: Highest ε4 frequency (46.1%)
- **Asian populations**: Higher ε2 frequency (19.2%)
- **European populations**: Moderate ε4 frequency (28.8%)

### Carrier Rates
- **ε4 carriers**: 27.3%
- **ε2 carriers**: 14.4%

## Adaptation for TCGA Cancer Studies

### Methodology
1. Use the same SNP coordinates (rs429358, rs7412)
2. Apply identical haplotype assignment logic
3. Stratify by cancer type instead of population
4. Test associations with cancer outcomes

### Potential Applications
- Compare APOE frequencies between cancer types
- Test APOE association with survival outcomes
- Analyze treatment response by APOE status
- Consider population ancestry in analyses

## Functions Overview

### Core Functions
- `verify_snp_coordinates()` - Define APOE SNP information
- `download_1000g_data()` - Download 1000 Genomes data
- `extract_snps_from_vcf()` - Extract target SNPs from VCF
- `build_haplotypes()` - Assign APOE isoforms from genotypes
- `compute_summaries()` - Calculate frequencies and carrier rates
- `perform_statistical_tests()` - Statistical analysis
- `create_visualizations()` - Generate plots
- `export_results()` - Save results to files

### Main Function
- `run_apoe_analysis()` - Complete workflow execution

## Data Requirements

### Input Format
- **VCF files**: Phased genotypes (0|1, 1|0, etc.)
- **SNP coordinates**: GRCh37/hg19 build
- **Sample metadata**: Population/cancer type information

### Quality Control
- Remove samples with missing genotypes
- Verify phasing format
- Check for non-canonical haplotypes

## Troubleshooting

### Common Issues
1. **Memory issues**: Use the simplified version (`run_apoe_simple.R`)
2. **Download failures**: Check internet connection and try again
3. **Package errors**: Install missing packages with `install.packages()`

### Data Verification
- Check that VCF files contain the target SNPs
- Verify sample metadata matches genotype data
- Ensure phased genotype format (0|1, 1|0)

## Clinical Relevance

### APOE and Disease Risk
- **ε4**: Major risk factor for Alzheimer's disease and cardiovascular disease
- **ε2**: Protective against Alzheimer's but increased cardiovascular risk
- **ε3**: Most common, considered "neutral"

### Population Health Implications
- African populations have highest ε4 frequency → higher Alzheimer's risk
- Asian populations have higher ε2 frequency → lower Alzheimer's risk
- Important for personalized medicine and risk assessment

## References

1. 1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. Nature, 526(7571), 68-74.
2. Corder, E. H., et al. (1993). Gene dose of apolipoprotein E type 4 allele and the risk of Alzheimer's disease in late onset families. Science, 261(5123), 921-923.
3. Mahley, R. W., & Rall, S. C. (2000). Apolipoprotein E: far more than a lipid transport protein. Annual review of genomics and human genetics, 1(1), 507-537.

## Support

For questions or issues:
1. Check the comprehensive documentation in `APOE_ANALYSIS_SUMMARY.md`
2. Review the code comments in `apoe_analysis_complete.R`
3. Verify data format and requirements

## License

This workflow is provided for research purposes. Please cite appropriately when using in publications.
