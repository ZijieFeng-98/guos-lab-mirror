# Secure Deployment Guide

## 🔒 Security Best Practices

### Never Commit Credentials to Git
- **NEVER** include API keys, tokens, or secrets in your code
- **NEVER** commit `.Rprofile` files with credentials
- **ALWAYS** use environment variables or secure configuration files

## 🛡️ How to Deploy Securely

### Option 1: Local .Rprofile (Recommended)
1. Create a `.Rprofile` file in your home directory (`~/.Rprofile`)
2. Add your ShinyApps.io credentials:
   ```r
   rsconnect::setAccountInfo(
     name   = 'your-username',
     token  = 'your-token',
     secret = 'your-secret'
   )
   ```
3. The `.Rprofile` file is automatically ignored by Git

### Option 2: Environment Variables
Set environment variables before running R:
```bash
export SHINYAPPS_NAME="your-username"
export SHINYAPPS_TOKEN="your-token"
export SHINYAPPS_SECRET="your-secret"
R
```

### Option 3: GitHub Secrets (for GitHub Actions)
1. Go to your GitHub repository
2. Settings > Secrets and variables > Actions
3. Add these secrets:
   - `SHINYAPPS_NAME`
   - `SHINYAPPS_TOKEN`
   - `SHINYAPPS_SECRET`

## 🚀 Deployment Commands

### Local Deployment
```r
# Your credentials are automatically loaded from .Rprofile
source("deploy_qpcr_app.R")
```

### Manual Deployment
```r
# Set credentials manually (not recommended for scripts)
rsconnect::setAccountInfo(
  name   = 'your-username',
  token  = 'your-token',
  secret = 'your-secret'
)

# Deploy
rsconnect::deployApp(appName = "qpcr-cleaner")
```

## ⚠️ Security Checklist

- [ ] No credentials in code files
- [ ] No credentials in commit history
- [ ] .Rprofile is in .gitignore
- [ ] Environment variables used for automation
- [ ] GitHub secrets configured (if using GitHub Actions)

## 🔍 Verify Security

Check that your credentials are not in any tracked files:
```bash
git grep -i "token\|secret\|password" -- "*.R" "*.md" "*.yml"
```

This should return no results if credentials are properly secured.
