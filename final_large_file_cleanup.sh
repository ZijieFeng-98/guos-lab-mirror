#!/bin/zsh
set -euo pipefail

echo "==> Final cleanup to remove remaining large files (>100MB)"

# Set environment variable to suppress git-filter-branch warnings
export FILTER_BRANCH_SQUELCH_WARNING=1

# Remove the remaining large files that exceed GitHub's 100MB limit
echo "==> Removing files >100MB from Git history"
git filter-branch --force --index-filter \
  'git rm --cached --ignore-unmatch \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/processed/GSE222520_combined_expression_matrix.rds" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/processed/GSE222520_combined_expression_matrix_test.csv" \
    "My_Research/Immunology_analysis/results/IMMUNOLOGY_ANALYSIS_CLEAN/DATA/paad_data.rds" \
    "My_Research/Immunology_analysis/results/IMMUNOLOGY_ANALYSIS_CLEAN/DATA/ov_data.rds" \
    "My_Research/Immunology_analysis/results/IMMUNOLOGY_ANALYSIS_CLEAN/DATA/gbm_data_new.rds"' \
  --prune-empty --tag-name-filter cat -- --all

# Clean up the repository
echo "==> Cleaning up repository"
git for-each-ref --format="%(refname)" refs/original/ | xargs -n 1 git update-ref -d
git reflog expire --expire=now --all
git gc --prune=now --aggressive

echo "==> Repository size after final cleanup:"
du -sh .git

echo "==> Checking for any remaining large files:"
find . -type f -size +100M -not -path "./.git/*" || echo "No files >100MB found"

echo "==> Ready to push (should work now)"
