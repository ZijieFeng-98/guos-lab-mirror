#!/bin/zsh
set -euo pipefail

echo "==> Aggressive cleanup of mirror repository to get under 2GB"

# Remove all large data files and archives from Git history
echo "==> Removing large data files from Git history"
git filter-branch --force --index-filter \
  'git rm --cached --ignore-unmatch \
    "My_Research/Immunology_analysis/data/raw/brca_data.rds" \
    "My_Research/Immunology_analysis/data/raw/gbm_data_new.rds" \
    "My_Research/Immunology_analysis/data/raw/ov_data.rds" \
    "My_Research/Immunology_analysis/data/raw/paad_data.rds" \
    "My_Research/Immunology_analysis/data/raw/archives/Mon_Aug_18_10_36_31_2025_0.tar.gz" \
    "My_Research/Immunology_analysis/data/raw/archives/Mon_Aug_18_10_36_31_2025_1.tar.gz" \
    "My_Research/Project_ApoE/APOE analysis/data/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
    "My_Research/Project_ApoE/APOE analysis/data/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSE222520_raw.rds"' \
  --prune-empty --tag-name-filter cat -- --all

# Remove all TCGA data files (they're very large)
echo "==> Removing TCGA data files"
git filter-branch --force --index-filter \
  'git rm --cached --ignore-unmatch \
    "My_Research/Immunology_analysis/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/" \
    "My_Research/Immunology_analysis/GDCdata/TCGA-GBM/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/" \
    "My_Research/Immunology_analysis/GDCdata/TCGA-OV/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/" \
    "My_Research/Immunology_analysis/GDCdata/TCGA-PAAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/"' \
  --prune-empty --tag-name-filter cat -- --all

# Clean up the repository
echo "==> Cleaning up repository"
git for-each-ref --format="%(refname)" refs/original/ | xargs -n 1 git update-ref -d
git reflog expire --expire=now --all
git gc --prune=now --aggressive

echo "==> Repository size after cleanup:"
du -sh .git

echo "==> Ready to push (if under 2GB)"
