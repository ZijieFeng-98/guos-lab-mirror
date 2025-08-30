# Quick Start Guide - qPCR Cleaner and Analyzer

## 🚀 Get Started in 5 Minutes

### Step 1: Install R and RStudio
1. Download and install R from [https://cran.r-project.org/](https://cran.r-project.org/)
2. Download and install RStudio from [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/)

### Step 2: Install Required Packages
Open R or RStudio and run:

```r
# Run the installation script
source("install_packages.R")
```

Or install manually:
```r
install.packages(c("shiny", "gtools", "readxl", "DT", "ggplot2", "dplyr", "tidyr", "shinyjs", "shinyWidgets"))
```

### Step 3: Run the App
```r
# Navigate to the app directory (if needed)
setwd("path/to/qPCR_Cleaner_Analyzer")

# Run the app
shiny::runApp()
```

### Step 4: Test with Sample Data
1. Upload the included `sample_data.csv` file
2. Select "3 replicates" 
3. Click "Run Analysis"
4. Explore the results!

## 📊 Sample Data Format

The included `sample_data.csv` contains:
- **10 samples** (rows)
- **3 genes** with **3 replicates each** (9 columns total)
- **Column structure**: 
  - Columns 1-3: Reference gene
  - Columns 4-6: Target gene 1  
  - Columns 7-9: Target gene 2

## 🎯 Basic Workflow

1. **Upload** your qPCR data file (CSV or Excel)
2. **Configure** analysis settings (replicates, outlier removal)
3. **Run** the analysis
4. **Review** results in the interactive tabs
5. **Download** complete results as CSV

## 🔧 Common Settings

- **3 replicates**: Standard qPCR setup
- **4 replicates**: Enable outlier removal for better quality control
- **Outlier threshold**: 2.0 SD (default), adjust if needed

## 📈 What You'll Get

- **Cleaned Ct values** with outliers removed
- **ΔCt calculations** with optimal replicate pairing
- **Expression values** (2^ΔCt)
- **Relative expression** compared to control
- **Statistical summaries** (mean, SD, CV)
- **Quality control metrics**
- **Interactive visualizations**

## 🆘 Need Help?

- Check the full [README.md](README.md) for detailed documentation
- Review the "Instructions" tab in the app
- Use the sample data to test functionality
- Ensure your data follows the required format

---

**Ready to analyze your qPCR data? Start with the sample data to see how it works!**
