# GitHub Setup Guide - qPCR Cleaner and Analyzer

This guide will help you set up your qPCR Cleaner and Analyzer project on GitHub and deploy it to ShinyApps.io.

## 🚀 Quick Setup

### 1. Automated Setup
Run the setup script to automatically configure your repository:

```r
source("setup_github.R")
```

This script will:
- Update placeholder text in files
- Initialize Git repository
- Create initial commit
- Provide step-by-step GitHub setup instructions

### 2. Manual Setup
If you prefer manual setup, follow the steps below.

## 📋 Prerequisites

- Git installed on your system
- GitHub account
- ShinyApps.io account (for deployment)
- R and RStudio

## 🔧 Step-by-Step Setup

### Step 1: Create GitHub Repository

1. Go to [GitHub](https://github.com/new)
2. Repository name: `qPCR_Cleaner_Analyzer`
3. Description: `qPCR Cleaner and Analyzer - A comprehensive Shiny app for qPCR data analysis`
4. Make it **Public**
5. **Don't** initialize with README (we already have one)
6. Click "Create repository"

### Step 2: Update Configuration Files

Update the following files with your information:

#### README.md
Replace `yourusername` with your actual GitHub username in:
- Badge URLs
- Installation instructions
- Support links

#### DESCRIPTION
Update the author information:
```r
Authors@R: 
    person("Your", "Name", 
           email = "your-email@example.com",
           role = c("aut", "cre"))
```

#### deploy_qpcr_app.R
Update the ShinyApps.io credentials:
```r
rsconnect::setAccountInfo(
  name   = 'your-shinyapps-username',
  token  = 'your-token',
  secret = 'your-secret'
)
```

### Step 3: Initialize Git Repository

```bash
# Initialize Git repository
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: qPCR Cleaner and Analyzer"

# Add remote origin
git remote add origin https://github.com/YOUR_USERNAME/qPCR_Cleaner_Analyzer.git

# Set main branch
git branch -M main

# Push to GitHub
git push -u origin main
```

### Step 4: Configure GitHub Repository

#### Add Repository Topics
Go to your repository and add these topics:
- `r`
- `shiny`
- `qpcr`
- `bioinformatics`
- `data-analysis`
- `statistics`
- `research`

#### Enable GitHub Pages (Optional)
1. Go to Settings > Pages
2. Source: Deploy from a branch
3. Branch: main, folder: / (root)
4. Save

#### Set up GitHub Actions (Optional)
1. Go to Settings > Secrets and variables > Actions
2. Add the following secrets:
   - `SHINYAPPS_NAME`: Your ShinyApps.io username
   - `SHINYAPPS_TOKEN`: Your ShinyApps.io token
   - `SHINYAPPS_SECRET`: Your ShinyApps.io secret

## 🚀 Deployment Options

### Option 1: Local Deployment
```r
# Run locally
shiny::runApp()

# Or use the package function
qPCR_Cleaner_Analyzer::run_app()
```

### Option 2: ShinyApps.io Deployment
```r
# Install rsconnect
install.packages("rsconnect")

# Configure your account (one-time setup)
rsconnect::setAccountInfo(
  name   = 'your-username',
  token  = 'your-token',
  secret = 'your-secret'
)

# Deploy using the deployment script
source("deploy_qpcr_app.R")
```

### Option 3: GitHub Actions (Automatic)
The repository includes GitHub Actions for automatic deployment. When you push to the main branch, it will automatically deploy to ShinyApps.io.

## 📦 Package Installation

### From GitHub
```r
# Install devtools if not already installed
if (!require(devtools)) install.packages("devtools")

# Install from GitHub
devtools::install_github("yourusername/qPCR_Cleaner_Analyzer")

# Run the app
qPCR_Cleaner_Analyzer::run_app()
```

### From Local Repository
```r
# Install from local directory
devtools::install(".")

# Run the app
qPCR_Cleaner_Analyzer::run_app()
```

## 🔍 Repository Structure

```
qPCR_Cleaner_Analyzer/
├── .github/
│   ├── workflows/
│   │   └── deploy.yml          # GitHub Actions workflow
│   └── ISSUE_TEMPLATE/
│       ├── bug_report.md       # Bug report template
│       └── feature_request.md  # Feature request template
├── R/
│   └── run_app.R              # Package function to run app
├── inst/
│   └── app/
│       ├── app.R              # Main Shiny application
│       └── sample_data.csv    # Sample data for testing
├── app.R                      # Main app file (for direct use)
├── DESCRIPTION                # R package description
├── NAMESPACE                  # R package namespace
├── LICENSE                    # MIT License
├── README.md                  # Main documentation
├── CONTRIBUTING.md            # Contributing guidelines
├── QUICK_START.md             # Quick start guide
├── GITHUB_SETUP.md            # This file
├── install_packages.R         # Package installation script
├── deploy_qpcr_app.R          # Deployment script
├── setup_github.R             # GitHub setup script
├── run_app.bat                # Windows launcher
├── sample_data.csv            # Sample data
└── .gitignore                 # Git ignore file
```

## 🛠️ Maintenance

### Updating the App
```bash
# Make your changes
git add .
git commit -m "Description of changes"
git push origin main
```

### Adding New Features
1. Create a feature branch: `git checkout -b feature/new-feature`
2. Make your changes
3. Test thoroughly
4. Commit: `git commit -m "Add new feature"`
5. Push: `git push origin feature/new-feature`
6. Create a Pull Request

### Version Updates
1. Update version in `DESCRIPTION`
2. Update `README.md` if needed
3. Create a release on GitHub
4. Tag the release

## 📞 Support

### For Users
- Create an issue using the bug report template
- Check the README.md for documentation
- Use the sample data for testing

### For Contributors
- Read CONTRIBUTING.md
- Follow the code style guidelines
- Test with sample data
- Update documentation

## 🎉 Success!

Once setup is complete, your qPCR Cleaner and Analyzer will be:
- ✅ Available on GitHub
- ✅ Installable as an R package
- ✅ Deployable to ShinyApps.io
- ✅ Ready for community contributions
- ✅ Professionally documented

Your app will be accessible at:
- **GitHub**: https://github.com/yourusername/qPCR_Cleaner_Analyzer
- **ShinyApps.io**: https://yourusername.shinyapps.io/qpcr-cleaner
- **Install**: `devtools::install_github("yourusername/qPCR_Cleaner_Analyzer")`

---

**Need help?** Check the troubleshooting section in README.md or create an issue on GitHub.
