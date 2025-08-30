# qPCR Cleaner and Analyzer

[![R](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://cran.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-1.7%2B-orange.svg)](https://shiny.rstudio.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/ZijieFeng-98/qPCR_Cleaner_Analyzer.svg)](https://github.com/ZijieFeng-98/qPCR_Cleaner_Analyzer/stargazers)
[![GitHub issues](https://img.shields.io/github/issues/ZijieFeng-98/qPCR_Cleaner_Analyzer.svg)](https://github.com/ZijieFeng-98/qPCR_Cleaner_Analyzer/issues)

A comprehensive Shiny web application for processing and analyzing qPCR (quantitative polymerase chain reaction) data with advanced statistical methods and interactive visualizations.

## 🚀 Features

### Core Analysis
- **Data Import**: Support for CSV and Excel files with or without headers
- **Flexible Replicate Handling**: Process 3 or 4 replicates per gene
- **Outlier Detection**: Automatic outlier identification and removal for 4-replicate data
- **ΔCt Calculation**: Optimal replicate pairing to minimize standard deviation
- **Expression Analysis**: Convert ΔCt to expression values (2^ΔCt)
- **Relative Expression**: Calculate relative expression compared to control samples
- **Statistical Summaries**: Mean, standard deviation, and coefficient of variation

### Quality Control
- **Data Validation**: Automatic detection of missing values and data quality issues
- **Outlier Analysis**: Visual and statistical outlier detection
- **Replicate Variability**: Coefficient of variation analysis across replicates
- **Quality Metrics**: Comprehensive quality control reporting

### Visualization
- **Interactive Tables**: Sortable and searchable data tables
- **Box Plots**: Distribution analysis of Ct values and expression data
- **Heatmaps**: Relative expression visualization across samples and genes
- **CV Plots**: Replicate variability analysis with quality thresholds

### User Experience
- **Intuitive Interface**: Tabbed layout with clear navigation
- **Progress Indicators**: Real-time analysis progress updates
- **Error Handling**: Comprehensive error messages and validation
- **Export Options**: Download complete analysis results as CSV

## 📋 Requirements

### R Packages
```r
# Core packages
install.packages(c("shiny", "gtools", "readxl"))

# Enhanced UI and visualization packages
install.packages(c("DT", "ggplot2", "dplyr", "tidyr", "shinyjs", "shinyWidgets"))
```

### System Requirements
- R version 3.6.0 or higher
- RStudio (recommended) or any R environment
- Internet connection for package installation

## 🛠️ Installation

### Option 1: Install from GitHub
```r
# Install devtools if not already installed
if (!require(devtools)) install.packages("devtools")

# Install from GitHub
devtools::install_github("ZijieFeng-98/qPCR_Cleaner_Analyzer")

# Run the app
qPCR_Cleaner_Analyzer::run_app()
```

### Option 2: Run Directly
```r
# Clone and run directly
git clone https://github.com/ZijieFeng-98/qPCR_Cleaner_Analyzer.git
cd qPCR_Cleaner_Analyzer
source("scripts/install_packages.R")
shiny::runApp()
```

### Option 3: Clone Repository
```bash
# Clone the repository
git clone https://github.com/ZijieFeng-98/qPCR_Cleaner_Analyzer.git
cd qPCR_Cleaner_Analyzer

# Install packages and run
Rscript scripts/install_packages.R
Rscript -e "shiny::runApp()"
```

### Option 4: Direct Download
1. Download the repository as ZIP
2. Extract to your desired location
3. Open R/RStudio in the extracted folder
4. Run the installation script:
   ```r
   source("scripts/install_packages.R")
   shiny::runApp()
   ```

### Option 5: Windows One-Click
- Double-click `scripts/run_app.bat` (Windows only)
- The script will automatically install packages and start the app

## 📊 Data Format Requirements

### File Format
- **Supported formats**: CSV (.csv) or Excel (.xlsx, .xls)
- **Headers**: Optional - can be enabled/disabled in the app
- **Missing values**: Can be marked as "Undetermined", "#VALUE!", "NA", or left blank

### Data Organization
- **Rows**: Individual samples/experiments
- **Columns**: Genes (organized in replicate blocks)
- **First gene**: Must be your reference/control gene
- **Replicate structure**: 
  - 3 replicates: Columns 1-3 = Gene1, Columns 4-6 = Gene2, etc.
  - 4 replicates: Columns 1-4 = Gene1, Columns 5-8 = Gene2, etc.

### Example Data Structure
```
Sample1: 20.5, 20.3, 20.7, 25.1, 25.3, 24.9, 30.2, 30.5, 30.1
Sample2: 21.1, 20.9, 21.3, 26.2, 26.0, 26.4, 31.1, 31.3, 30.9
```

Where:
- Columns 1-3: Reference gene (3 replicates)
- Columns 4-6: Target gene 1 (3 replicates)  
- Columns 7-9: Target gene 2 (3 replicates)

## 🎯 Usage Guide

### Step 1: Upload Data
1. Click "Choose File" in the Data Upload section
2. Select your CSV or Excel file
3. Check "File has column headers" if applicable
4. Verify data preview in the "Data Preview" tab

### Step 2: Configure Analysis
1. Select number of replicates per gene (3 or 4)
2. For 4 replicates, enable outlier removal if desired
3. Adjust outlier detection threshold (default: 2.0 SD)

### Step 3: Run Analysis
1. Click "🚀 Run Analysis" button
2. Monitor progress in the progress bar
3. Review results in the "Analysis Results" tab

### Step 4: Review Results
1. **Summary Statistics**: View mean, SD, and CV for each gene
2. **Detailed Results**: Examine cleaned Ct values, ΔCt, expression, and relative expression
3. **Visualizations**: Explore data distributions and patterns
4. **Quality Control**: Check data quality metrics and outlier detection

### Step 5: Export Results
1. Click "📥 Download Results" to save complete analysis as CSV
2. File includes all analysis sections with clear headers

## 📈 Analysis Methods

### Outlier Detection
- **Method**: Standard deviation-based detection
- **Threshold**: Configurable (default: 2.0 SD from mean)
- **Application**: Only for 4-replicate data
- **Process**: Removes most extreme outlier per gene block

### ΔCt Calculation
- **Formula**: ΔCt = Ct(reference) - Ct(target)
- **Optimization**: Permutation-based pairing to minimize SD
- **Fallback**: Simple pairing if optimization fails

### Expression Analysis
- **Formula**: Expression = 2^ΔCt
- **Purpose**: Convert ΔCt to relative expression values
- **Interpretation**: Higher values indicate higher expression

### Relative Expression
- **Formula**: Relative = Expression(treatment) / Expression(control)
- **Control**: First sample (row 1) used as reference
- **Optimization**: Permutation-based pairing for minimal variance

### Quality Metrics
- **Coefficient of Variation**: CV% = (SD/Mean) × 100
- **Quality Threshold**: CV% < 20% recommended
- **Missing Data**: Percentage of missing values reported

## 🔧 Advanced Features

### Customizable Outlier Detection
- Adjustable threshold (1.0-5.0 SD)
- Visual outlier identification
- Detailed outlier reporting

### Interactive Visualizations
- **Box Plots**: Distribution analysis with individual data points
- **Heatmaps**: Color-coded relative expression patterns
- **CV Plots**: Quality control with threshold indicators

### Comprehensive Export
- All analysis sections included
- Formatted with clear headers
- Compatible with Excel and statistical software

## 🚀 Deployment

### Local Deployment
```r
# Run locally
shiny::runApp()

# Run with custom port
shiny::runApp(port = 8080)
```

### ShinyApps.io Deployment
```r
# Install rsconnect
install.packages("rsconnect")

# Configure your account (one-time setup)
rsconnect::setAccountInfo(name='yourname', token='yourtoken', secret='yoursecret')

# Deploy
source("scripts/deploy_qpcr_app.R")
```

### GitHub Actions (Automatic)
The repository includes GitHub Actions for automatic deployment to ShinyApps.io. Set up the following secrets in your repository:
- `SHINYAPPS_NAME`: Your ShinyApps.io username
- `SHINYAPPS_TOKEN`: Your ShinyApps.io token
- `SHINYAPPS_SECRET`: Your ShinyApps.io secret

## 🐛 Troubleshooting

### Common Issues

**"Uploaded file has no usable data after cleaning"**
- Check file format (CSV or Excel)
- Verify data contains numeric values
- Ensure file is not empty or corrupted

**"Analysis error"**
- Verify correct number of replicates selected
- Check that first gene is reference gene
- Ensure sufficient data for analysis

**Missing packages**
```r
# Install all required packages
install.packages(c("shiny", "gtools", "readxl", "DT", "ggplot2", "dplyr", "tidyr", "shinyjs", "shinyWidgets"))
```

### Performance Tips
- For large datasets (>1000 samples), consider processing in batches
- Close other applications to free up memory
- Use CSV format for faster loading

## 📝 Output Files

The downloaded CSV file contains the following sections:

1. **Cleaned Ct Data**: Original data with outliers removed
2. **ΔCt Values**: Calculated ΔCt values with optimal pairing
3. **Expression Values**: 2^ΔCt conversion
4. **Expression Means**: Mean expression per gene per sample
5. **Expression SDs**: Standard deviation of expression
6. **Relative Expression**: Expression relative to control
7. **Relative Expression Means**: Mean relative expression
8. **Relative Expression SDs**: Standard deviation of relative expression

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### How to Contribute
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Reporting Issues
- Use the [bug report template](.github/ISSUE_TEMPLATE/bug_report.md)
- Include sample data and error messages
- Specify your environment details

### Feature Requests
- Use the [feature request template](.github/ISSUE_TEMPLATE/feature_request.md)
- Describe the use case and benefits
- Consider implementation complexity

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Built with R Shiny framework
- Statistical methods based on standard qPCR analysis protocols
- Visualization powered by ggplot2
- Interactive tables by DataTables
- Community contributors and feedback

## 📞 Support

- 📧 **Email**: [fengzijiehh@outlook.com]
- 🐛 **Issues**: [GitHub Issues](https://github.com/ZijieFeng-98/qPCR_Cleaner_Analyzer/issues)
- 📖 **Documentation**: [Wiki](https://github.com/ZijieFeng-98/qPCR_Cleaner_Analyzer/wiki)
- 💬 **Discussions**: [GitHub Discussions](https://github.com/ZijieFeng-98/qPCR_Cleaner_Analyzer/discussions)

---

**For questions or support, please refer to the troubleshooting section or create an issue in the repository.**
