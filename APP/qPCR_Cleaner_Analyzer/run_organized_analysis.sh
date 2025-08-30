#!/bin/bash

# Immunology Analysis Project - Organized Analysis Runner
# This script helps you run analysis from the organized project structure

echo "🧬 Immunology Analysis Project - Organized Analysis Runner"
echo "=========================================================="
echo ""

# Check if we're in the right directory
if [ ! -f "Immunology_analysis.Rproj" ]; then
    echo "❌ Error: Please run this script from the Immunology_analysis project root directory"
    exit 1
fi

# Function to display menu
show_menu() {
    echo ""
    echo "📋 Available Analysis Options:"
    echo "1. Install/Check Dependencies"
    echo "2. Run GBM Analysis"
    echo "3. Run BRCA Analysis"
    echo "4. Run OV Analysis"
    echo "5. Run PAAD Analysis"
    echo "6. Run All Cancer Types"
    echo "7. Check GDC Server Status"
    echo "8. View Project Structure"
    echo "9. Open Results Directory"
    echo "0. Exit"
    echo ""
}

# Function to run analysis
run_analysis() {
    local cancer_type=$1
    local script_name=$2
    
    echo "🔬 Running $cancer_type analysis..."
    echo "Script: $script_name"
    echo ""
    
    if [ -f "scripts/analysis/$script_name" ]; then
        Rscript "scripts/analysis/$script_name"
        echo ""
        echo "✅ $cancer_type analysis completed!"
        echo "📊 Results saved to results/ directory"
    else
        echo "❌ Error: Script not found: scripts/analysis/$script_name"
    fi
}

# Main menu loop
while true; do
    show_menu
    read -p "Enter your choice (0-9): " choice
    
    case $choice in
        1)
            echo "📦 Installing/Checking Dependencies..."
            Rscript scripts/utilities/install_packages.R
            ;;
        2)
            run_analysis "GBM" "TCGA_GBM_Correlation_Analysis_Expanded_Immune.R"
            ;;
        3)
            run_analysis "BRCA" "TCGA_BRCA_Correlation_Analysis.R"
            ;;
        4)
            run_analysis "OV" "TCGA_OV_Correlation_Analysis.R"
            ;;
        5)
            run_analysis "PAAD" "TCGA_PAAD_Correlation_Analysis.R"
            ;;
        6)
            echo "🔬 Running All Cancer Type Analyses..."
            echo "This may take a while..."
            echo ""
            
            run_analysis "GBM" "TCGA_GBM_Correlation_Analysis_Expanded_Immune.R"
            run_analysis "BRCA" "TCGA_BRCA_Correlation_Analysis.R"
            run_analysis "OV" "TCGA_OV_Correlation_Analysis.R"
            run_analysis "PAAD" "TCGA_PAAD_Correlation_Analysis.R"
            
            echo "🎉 All analyses completed!"
            ;;
        7)
            echo "🌐 Checking GDC Server Status..."
            Rscript scripts/utilities/check_gdc_status.R
            ;;
        8)
            echo "📁 Project Structure:"
            echo ""
            echo "scripts/"
            echo "├── analysis/          # Main analysis scripts"
            echo "├── data_download/     # Data download utilities"
            echo "└── utilities/         # Setup and utility scripts"
            echo ""
            echo "data/"
            echo "├── raw/              # Raw TCGA data files"
            echo "└── processed/        # Processed data"
            echo ""
            echo "results/"
            echo "├── figures/          # Generated plots"
            echo "└── tables/           # Statistical results"
            echo ""
            echo "docs/                 # Documentation"
            echo ""
            echo "📖 For detailed structure, see: docs/PROJECT_STRUCTURE.md"
            ;;
        9)
            echo "📂 Opening Results Directory..."
            if [[ "$OSTYPE" == "darwin"* ]]; then
                open results/
            elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
                xdg-open results/
            else
                echo "Please manually open the results/ directory"
            fi
            ;;
        0)
            echo "👋 Goodbye!"
            exit 0
            ;;
        *)
            echo "❌ Invalid choice. Please enter a number between 0-9."
            ;;
    esac
    
    echo ""
    read -p "Press Enter to continue..."
done

