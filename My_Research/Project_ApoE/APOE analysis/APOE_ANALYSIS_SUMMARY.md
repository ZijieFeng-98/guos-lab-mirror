# APOE Isoform Analysis - Complete Summary

## Overview
This document summarizes the complete APOE (Apolipoprotein E) isoform analysis workflow, including the underlying logic, methodology, and results from the 1000 Genomes Project data analysis.

## APOE Biology and Clinical Relevance

### APOE Gene and Protein
- **Gene**: APOE (Apolipoprotein E)
- **Location**: Chromosome 19 (GRCh37/hg19)
- **Function**: Lipid transport and metabolism
- **Clinical Significance**: Major genetic risk factor for Alzheimer's disease and cardiovascular disease

### APOE Isoforms
APOE has three major isoforms determined by two SNPs:

1. **ε2 (E2)**: Protective against Alzheimer's disease, but increased cardiovascular risk
2. **ε3 (E3)**: Most common, considered "neutral" 
3. **ε4 (E4)**: Major risk factor for Alzheimer's disease and cardiovascular disease

## Genetic Basis of APOE Isoforms

### Key SNPs
| SNP ID | Chromosome | Position (GRCh37) | Reference Allele | Alternate Allele |
|--------|------------|-------------------|------------------|------------------|
| rs429358 | 19 | 45411941 | T | C |
| rs7412 | 19 | 45412079 | C | T |

### Haplotype Assignment Logic
From phased genotypes (e.g., 0|1, 1|0), we construct haplotypes:

**Haplotype Construction:**
- Haplotype 1: rs429358_allele1 + rs7412_allele1
- Haplotype 2: rs429358_allele2 + rs7412_allele2

**APOE Isoform Assignment:**
```
rs429358 | rs7412 | Haplotype | APOE Isoform
---------|--------|-----------|-------------
T        | T      | T-T       | ε2
T        | C      | T-C       | ε3  
C        | C      | C-C       | ε4
C        | T      | C-T       | ε1 (non-canonical)
```

**Genotype Assignment:**
- Individual genotypes are combinations of two haplotypes
- Sorted alphabetically (e.g., ε2/ε4 not ε4/ε2)
- Examples: ε2/ε2, ε2/ε3, ε2/ε4, ε3/ε3, ε3/ε4, ε4/ε4

**Carrier Status:**
- ε4 carrier: Has at least one ε4 allele
- ε2 carrier: Has at least one ε2 allele

## Analysis Workflow

### 1. Data Acquisition
- **Source**: 1000 Genomes Project Phase 3
- **Data**: Chromosome 19 VCF file (phased genotypes)
- **Samples**: 2,504 individuals from 5 superpopulations
- **Format**: VCF (Variant Call Format) with phased genotypes

### 2. Data Processing Steps
1. **SNP Extraction**: Extract rs429358 and rs7412 from VCF
2. **Genotype Extraction**: Extract phased genotypes (0|1, 1|0, etc.)
3. **Haplotype Construction**: Build haplotypes from phased genotypes
4. **APOE Assignment**: Assign ε2, ε3, ε4, ε1 based on haplotype logic
5. **Genotype Creation**: Combine haplotypes into individual genotypes
6. **Carrier Status**: Determine ε4 and ε2 carrier status

### 3. Quality Control
- Remove samples with missing genotypes
- Verify phasing (should be 0|1, 1|0 format)
- Check for non-canonical haplotypes (ε1)

## Results from 1000 Genomes Analysis

### Overall Frequencies

#### Haplotype Frequencies
| Haplotype | Count | Frequency | Percentage |
|-----------|-------|-----------|------------|
| ε3 | 3,879 | 0.775 | 77.5% |
| ε4 | 753 | 0.150 | 15.0% |
| ε2 | 375 | 0.075 | 7.5% |
| ε1 | 1 | 0.0002 | 0.02% |

#### Genotype Frequencies
| Genotype | Count | Frequency | Percentage |
|----------|-------|-----------|------------|
| ε3/ε3 | 1,525 | 0.609 | 60.9% |
| ε3/ε4 | 547 | 0.218 | 21.8% |
| ε2/ε3 | 282 | 0.113 | 11.3% |
| ε4/ε4 | 70 | 0.028 | 2.8% |
| ε2/ε4 | 65 | 0.026 | 2.6% |
| ε2/ε2 | 14 | 0.006 | 0.6% |
| ε1/ε4 | 1 | 0.0004 | 0.04% |

#### Carrier Rates
- **ε4 carriers**: 27.3% (683/2,504 individuals)
- **ε2 carriers**: 14.4% (361/2,504 individuals)

### Population Differences

#### ε4 Carrier Rates by Population
| Population | Sample Size | ε4 Carrier Rate | Percentage |
|------------|-------------|-----------------|------------|
| AFR (African) | 661 | 0.461 | 46.1% |
| EUR (European) | 503 | 0.288 | 28.8% |
| AMR (American) | 347 | 0.202 | 20.2% |
| EAS (East Asian) | 504 | 0.165 | 16.5% |
| SAS (South Asian) | 489 | 0.164 | 16.4% |

#### ε2 Carrier Rates by Population
| Population | Sample Size | ε2 Carrier Rate | Percentage |
|------------|-------------|-----------------|------------|
| EAS (East Asian) | 504 | 0.192 | 19.2% |
| AFR (African) | 661 | 0.197 | 19.7% |
| EUR (European) | 503 | 0.123 | 12.3% |
| AMR (American) | 347 | 0.092 | 9.2% |
| SAS (South Asian) | 489 | 0.082 | 8.2% |

### Statistical Tests

#### Chi-Square Tests
- **ε4 carrier by population**: χ² = 156.8, p < 0.001
- **ε2 carrier by population**: χ² = 45.2, p < 0.001

**Interpretation**: Highly significant differences in APOE allele frequencies across populations.

## Clinical Implications

### Alzheimer's Disease Risk
- **ε4 carriers**: 3-15x increased risk depending on number of ε4 alleles
- **ε2 carriers**: 0.5-0.7x reduced risk (protective)
- **ε3/ε3**: Reference (baseline risk)

### Cardiovascular Disease Risk
- **ε4 carriers**: Increased cardiovascular disease risk
- **ε2 carriers**: Increased cardiovascular disease risk (paradoxical)
- **ε3/ε3**: Baseline cardiovascular risk

### Population Health Implications
- **African populations**: Highest ε4 frequency → higher Alzheimer's risk
- **Asian populations**: Higher ε2 frequency → lower Alzheimer's risk but higher cardiovascular risk
- **European populations**: Moderate ε4 frequency → moderate Alzheimer's risk

## Adaptation for TCGA Cancer Analysis

### Methodology Adaptation
1. **Same SNP coordinates**: Use rs429358 and rs7412
2. **Same haplotype logic**: Apply identical APOE assignment rules
3. **Stratification**: Group by cancer type instead of population
4. **Controls**: Compare cancer cases vs. healthy controls
5. **Outcomes**: Test associations with cancer progression, survival, treatment response

### Potential Cancer Associations
- **ε4**: May influence cancer risk, progression, or treatment response
- **ε2**: May have protective or risk effects depending on cancer type
- **Population differences**: Consider ancestry in cancer studies

### Statistical Approaches for TCGA
1. **Case-control studies**: Compare APOE frequencies between cancer and controls
2. **Survival analysis**: Test APOE association with overall survival, progression-free survival
3. **Treatment response**: Analyze APOE association with chemotherapy/immunotherapy response
4. **Multi-variable models**: Adjust for age, sex, stage, treatment, ancestry

## Technical Implementation

### Required R Packages
```r
library(vcfR)      # VCF file reading
library(dplyr)     # Data manipulation
library(tidyr)     # Data reshaping
library(ggplot2)   # Visualization
library(readr)     # File reading
library(stringr)   # String manipulation
library(purrr)     # Functional programming
library(forcats)   # Factor manipulation
```

### Key Functions
1. `verify_snp_coordinates()`: Define APOE SNP information
2. `extract_snps_from_vcf()`: Extract target SNPs from VCF
3. `build_haplotypes()`: Construct haplotypes and assign APOE isoforms
4. `compute_summaries()`: Calculate frequencies and carrier rates
5. `perform_statistical_tests()`: Chi-square tests and odds ratios

### Data Requirements
- **VCF format**: Phased genotypes (0|1, 1|0, etc.)
- **SNP coordinates**: GRCh37/hg19 build
- **Sample metadata**: Population/cancer type information
- **Quality control**: Missing genotype handling

## File Outputs

### Generated Files
1. `apoe_genotypes.csv`: Individual-level genotype data
2. `haplotype_frequencies.csv`: Overall haplotype frequencies
3. `genotype_frequencies.csv`: Overall genotype frequencies
4. `superpopulation_summary.csv`: Population-stratified results
5. `statistical_tests.txt`: Statistical test results
6. `*.png`: Visualization plots

### Key Variables
- `hap1_apoe`, `hap2_apoe`: Individual haplotypes
- `genotype`: Combined genotype (e.g., "ε3/ε4")
- `carrier_e4`, `carrier_e2`: Binary carrier status
- `superpop`: Population group

## Conclusion

This APOE analysis workflow provides a robust, reproducible method for analyzing APOE isoform frequencies in genomic datasets. The methodology can be directly adapted for TCGA cancer studies by:

1. Using the same SNP coordinates and haplotype logic
2. Stratifying by cancer type instead of population
3. Testing associations with cancer outcomes
4. Considering population ancestry in analyses

The results demonstrate significant population differences in APOE allele frequencies, with important implications for disease risk across different ethnic groups. This framework provides a solid foundation for future cancer genomics studies involving APOE.

## References

1. 1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. Nature, 526(7571), 68-74.
2. Corder, E. H., et al. (1993). Gene dose of apolipoprotein E type 4 allele and the risk of Alzheimer's disease in late onset families. Science, 261(5123), 921-923.
3. Mahley, R. W., & Rall, S. C. (2000). Apolipoprotein E: far more than a lipid transport protein. Annual review of genomics and human genetics, 1(1), 507-537.
