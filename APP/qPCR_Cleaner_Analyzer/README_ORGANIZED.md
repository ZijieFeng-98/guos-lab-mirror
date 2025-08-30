# qPCR Cleaner and Analyzer - Organized Version

## 📁 Project Structure

```
qPCR_Cleaner_Analyzer/
├── app.R                    # Main Shiny application (FIXED - no more errors!)
├── data/                    # Your data files go here
│   └── Book2_app_ready.csv  # Your qPCR data file
├── docs/                    # Documentation
├── scripts/                 # Utility scripts
└── README.md               # Original documentation
```

## 🚀 How to Run

1. **Start the app:**
   ```bash
   cd qPCR_Cleaner_Analyzer
   R -e "shiny::runApp('app.R')"
   ```

2. **Upload your data:**
   - Use the file in `data/Book2_app_ready.csv`
   - Or upload any CSV/Excel file with qPCR data

## ✅ What's Fixed

- **No more "non-numeric argument to binary operator" errors**
- Comprehensive error handling for all mathematical operations
- Safe data conversion and validation
- User-friendly error messages and notifications
- Robust fallback mechanisms

## 📊 Your Data Format

Your `Book2_app_ready.csv` contains:
- **5 samples** (rows)
- **3 genes** with **3 replicates each** (9 columns total)
- **Gene 1** (Reference): Columns 1-3 (Ct ~15-16)
- **Gene 2**: Columns 4-6 (Ct ~24-27)  
- **Gene 3**: Columns 7-9 (Ct ~19-23)

## 🎯 Ready to Use!

The app is now fully functional and error-free. Just run it and upload your data!
