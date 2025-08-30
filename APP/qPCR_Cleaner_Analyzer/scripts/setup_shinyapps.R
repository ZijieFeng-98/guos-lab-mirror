#!/usr/bin/env Rscript
# ShinyApps.io Setup Script
# This script helps you configure your ShinyApps.io credentials

cat("Setting up ShinyApps.io credentials...\n\n")

# Instructions
cat("To deploy to ShinyApps.io, you need to:\n\n")

cat("1. Create a ShinyApps.io account:\n")
cat("   - Go to https://www.shinyapps.io/\n")
cat("   - Sign up for a free account\n\n")

cat("2. Get your account token:\n")
cat("   - Log in to ShinyApps.io\n")
cat("   - Go to Account > Tokens\n")
cat("   - Click 'Show' to reveal your token\n")
cat("   - Copy the token information\n\n")

cat("3. Configure your credentials:\n")
cat("   Run the following command in R:\n\n")

cat("rsconnect::setAccountInfo(\n")
cat("  name   = 'YOUR_SHINYAPPS_USERNAME',\n")
cat("  token  = 'YOUR_TOKEN',\n")
cat("  secret = 'YOUR_SECRET'\n")
cat(")\n\n")

cat("4. Deploy your app:\n")
cat("   source('deploy_qpcr_app.R')\n\n")

cat("Example:\n")
cat("rsconnect::setAccountInfo(\n")
cat("  name   = 'your-username',\n")
cat("  token  = 'your-token',\n")
cat("  secret = 'your-secret'\n")
cat(")\n\n")

cat("⚠️  Note: Keep your token and secret secure and don't share them publicly!\n")
