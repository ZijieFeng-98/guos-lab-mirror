#!/bin/bash

# TCGA GBM Correlation Analysis Runner Script
# This script sets up the environment and runs the analysis

echo "=========================================="
echo "TCGA GBM Correlation Analysis"
echo "Expanded Immune Marker List"
echo "=========================================="
echo ""

# Check if R is installed
if ! command -v R &> /dev/null; then
    echo "Error: R is not installed or not in PATH"
    echo "Please install R from https://cran.r-project.org/"
    exit 1
fi

echo "✓ R found: $(R --version | head -n 1)"
echo ""

# Create Processed_Data directory if it doesn't exist
mkdir -p "Processed_Data"

echo "Step 1: Installing required packages..."
echo "======================================"
Rscript install_packages.R

echo ""
echo "Step 2: Running correlation analysis..."
echo "====================================="
echo "This may take several minutes depending on your internet connection and system performance."
echo "The script will download TCGA-GBM data (~2-5 GB) and perform correlation analysis."
echo ""

# Run the main analysis script
Rscript TCGA_GBM_Correlation_Analysis_Expanded_Immune.R

echo ""
echo "=========================================="
echo "Analysis Complete!"
echo "=========================================="
echo ""
echo "Output files generated:"
echo "- Processed_Data/Spearman_Correlation_Heatmap_All_Immune.png"
echo "- Processed_Data/Spearman_Correlation_Matrix_All_Immune.csv"
echo ""
echo "You can find the results in the Processed_Data/ directory."
