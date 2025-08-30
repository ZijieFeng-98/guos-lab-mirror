# =============================================================================
# APOE ANALYSIS RUNNER SCRIPT
# =============================================================================
# 
# This script runs the complete APOE isoform analysis workflow.
# It will download data, process genotypes, assign APOE isoforms,
# and generate comprehensive results and visualizations.
#
# Usage: Rscript run_apoe_analysis.R
# =============================================================================

cat("=== APOE ISOFORM ANALYSIS WORKFLOW ===\n")
cat("Starting comprehensive APOE analysis...\n\n")

# Source the main analysis script
source("apoe_analysis_complete.R")

cat("=== ANALYSIS COMPLETE ===\n")
cat("Check the 'output/' directory for all results.\n")
cat("Check 'APOE_ANALYSIS_SUMMARY.md' for detailed documentation.\n")
