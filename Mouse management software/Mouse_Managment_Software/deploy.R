# Deployment script for Mouse Management Software
# Run this to deploy to shinyapps.io

# Install rsconnect if not already installed
if (!require(rsconnect, quietly = TRUE)) {
  install.packages("rsconnect")
}

# Set account info (replace with your actual credentials)
rsconnect::setAccountInfo(
  name = 'zijiefeng', 
  token = '73D52A6F9E64A708252DDAC86A116A09', 
  secret = 'ICyttOI6F9bPdLpS3iaOp208zJADDN9H7J+yS19F'
)

# Deploy the app
rsconnect::deployApp(
  appName = "mouse-management-software",
  appTitle = "Zijie's Great Mouse Management Software",
  appFiles = c(
    "app.R",
    "README.md",
    "QUICK_START.md",
    "install_packages.R",
    "sample_data.R",
    "mLIMS_Template.csv",
    "mLIMS_Cage_Template.csv"
  ),
  forceUpdate = TRUE
)

cat("✅ App deployed successfully!\n")
cat("Visit: https://zijiefeng.shinyapps.io/mouse-management-software/\n")
