# 🚀 Deployment Guide for Mouse Management Software

## Quick Deploy to shinyapps.io

### 1. One-Click Deployment
```r
# Run the deployment script
source("deploy.R")
```

### 2. Manual Deployment
```r
# Install rsconnect
install.packages("rsconnect")

# Set your account info
rsconnect::setAccountInfo(
  name = 'zijiefeng', 
  token = '73D52A6F9E64A708252DDAC86A116A09', 
  secret = 'ICyttOI6F9bPdLpS3iaOp208zJADDN9H7J+yS19F'
)

# Deploy
rsconnect::deployApp()
```

## 🔧 Production Setup (Optional)

### Database Configuration
To enable persistent storage, set these environment variables in your shinyapps.io dashboard:

**Required for PostgreSQL:**
- `DB_HOST` - Your database host (e.g., `your-db.region.rds.amazonaws.com`)
- `DB_NAME` - Database name (e.g., `mlims`)
- `DB_USER` - Database username (e.g., `mlims_app`)
- `DB_PASSWORD` - Database password
- `DB_PORT` - Database port (default: `5432`)
- `DB_SSLMODE` - SSL mode (default: `require`)

**Optional for S3 Backups:**
- `S3_BUCKET` - Your S3 bucket name
- `S3_REGION` - AWS region (e.g., `us-east-1`)
- `AWS_ACCESS_KEY_ID` - AWS access key
- `AWS_SECRET_ACCESS_KEY` - AWS secret key

### Setting Environment Variables in shinyapps.io

1. Go to your app dashboard: https://www.shinyapps.io/
2. Click on your app
3. Go to **Settings** tab
4. Scroll to **Environment Variables**
5. Add each variable with its value
6. Click **Save**

### Database Setup

**Option 1: Supabase (Recommended for small labs)**
1. Create account at https://supabase.com
2. Create new project
3. Get connection details from Settings → Database
4. Set environment variables

**Option 2: AWS RDS**
1. Create PostgreSQL instance in AWS RDS
2. Configure security groups to allow connections
3. Get connection details from RDS console
4. Set environment variables

**Option 3: Neon**
1. Create account at https://neon.tech
2. Create new project
3. Get connection string from dashboard
4. Parse and set environment variables

### S3 Backup Setup (Optional)

1. Create S3 bucket in AWS
2. Enable versioning on the bucket
3. Create IAM user with S3 access
4. Set environment variables

## 🔒 Security Best Practices

### Database Security
- ✅ Use TLS/SSL connections (`DB_SSLMODE=require`)
- ✅ Use dedicated database user (not admin)
- ✅ Restrict database access by IP if possible
- ✅ Rotate passwords regularly

### AWS Security
- ✅ Use IAM roles instead of access keys when possible
- ✅ Enable S3 bucket encryption
- ✅ Enable S3 bucket versioning
- ✅ Set up CloudTrail for audit logs

### App Security
- ✅ Set app to "Private" in shinyapps.io
- ✅ Invite specific users only
- ✅ Use strong passwords
- ✅ Enable 2FA on your account

## 📊 Monitoring & Maintenance

### App Monitoring
- Monitor app usage in shinyapps.io dashboard
- Check logs for errors
- Monitor database connections

### Database Maintenance
- Set up automated backups
- Monitor database performance
- Clean up old audit logs periodically

### Backup Strategy
- Daily automated backups to S3
- Weekly manual backups
- Test restore procedures monthly

## 🚨 Troubleshooting

### Common Issues

**App won't start:**
- Check package dependencies
- Verify environment variables
- Check app logs in shinyapps.io

**Database connection fails:**
- Verify credentials
- Check network connectivity
- Ensure database is running

**Backup fails:**
- Check S3 credentials
- Verify bucket permissions
- Check network connectivity

### Getting Help
- Check app logs in shinyapps.io dashboard
- Review R console output
- Check database connection status in Settings tab

## 🔄 Updates & Redeployment

### Updating the App
```r
# Make your changes to app.R
# Then redeploy
rsconnect::deployApp()
```

### Database Migrations
If you add new columns, add this to your app:
```r
# Add new column safely
DBI::dbExecute(pool, "ALTER TABLE mice ADD COLUMN IF NOT EXISTS new_field TEXT;")
```

### Backup Before Updates
Always create a backup before major updates:
1. Go to Settings tab
2. Click "☁️ Backup Now"
3. Deploy your changes
4. Test the app

## 📈 Scaling Considerations

### For Small Labs (< 1000 mice)
- Basic shinyapps.io plan is sufficient
- Supabase free tier works well
- S3 backups optional

### For Medium Labs (1000-10000 mice)
- Consider shinyapps.io paid plan
- Use dedicated database instance
- Enable S3 backups
- Monitor performance

### For Large Labs (> 10000 mice)
- Consider Posit Connect or ShinyProxy
- Use dedicated database cluster
- Implement data archiving
- Set up monitoring and alerting

## 🎯 Next Steps

1. **Deploy the app** using the deployment script
2. **Configure database** (optional but recommended)
3. **Set up S3 backups** (optional)
4. **Invite team members** to the app
5. **Train users** on the new features
6. **Monitor usage** and gather feedback

Your app is now production-ready with enterprise-grade features! 🎉
