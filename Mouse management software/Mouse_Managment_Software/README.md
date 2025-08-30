# Zijie's Great Mouse Management Software

A comprehensive R Shiny application for managing mouse colonies with advanced UI/UX features, keyboard shortcuts, comprehensive filtering, gene catalog management, and accessibility improvements. Built for efficient mouse colony management with modern interface design.

## 🎯 Features Implemented

### ✅ 1. Clear Click Targets
- **Large, highlighted buttons** with hover effects
- **Card-based interface** for cages that highlight on hover
- **Clear labeling** for all interactive elements (View Details, Edit, Download)
- **Icons with tooltips** for better user guidance

### ✅ 2. Cursor Feedback System
| Action | Cursor Type | CSS Class |
|--------|-------------|-----------|
| Clickable elements | Pointer (🖱️) | `cursor: pointer` |
| Non-clickable text | Default | `cursor: default` |
| Drag & drop areas | Move (↕️) | `cursor: move` |
| Loading states | Wait (⏳) | `cursor: wait` |
| Invalid actions | Not allowed (🚫) | `cursor: not-allowed` |

### ✅ 3. Hover Highlights
- **Cage cards** turn blue with shadow on hover
- **Buttons** change color and lift slightly
- **Mouse IDs** become underlined and change color
- **Smooth transitions** for all hover effects

### ✅ 4. Keyboard Navigation
- **Tab navigation** through all interactive elements
- **Enter/Space** to activate buttons and cards
- **Arrow keys** for dropdown selections
- **Accessible focus indicators**

### ✅ 5. Enhanced Keyboard Shortcuts
| Action | Shortcut | Description |
|--------|----------|-------------|
| New Cage | `Alt + N` | Opens New Cage entry form |
| Download All Data | `Alt + D` | Triggers full Excel/CSV download |
| Print All | `Alt + P` | Opens print preview for all cages |
| Search/Filter Cages | `/` | Focuses on search bar (when not editing) |
| Add Mouse to Cage | `Alt + M` | Opens Add Mouse modal in cage view |
| Export Current Cage | `Alt + E` | Download Excel for current cage |
| Print Current Cage | `Alt + Shift + P` | Opens print for current cage only |
| Edit Mouse Record | `Alt + Shift + E` | Opens editable fields for selected |
| Upload Data | `Alt + U` | Opens data upload modal |
| Back to Home | `Esc` or `Alt + H` | Goes back to homepage |

**Smart Keyboard Handling**: Shortcuts are disabled when typing in input fields, allowing users to type symbols like `flox/+` without interference.

### ✅ 6. Comprehensive Cage Filtering System
- **Home Dashboard Filtering**: Advanced filters for cage overview table
- **Cage View List Filtering**: Identical filtering capabilities in cage list view
- **Filter Options**:
  - **Gene/Genotype**: Filter by specific genes and genotype symbols
  - **Sex**: Filter by M, F, Mixed, or Any
  - **Status**: Filter by Breeding, Holding, Experiment, Active, Retired
  - **Has Pups**: Filter cages with or without pups
  - **Genotyping Status**: Filter by genotyping done/confirmed status
  - **Setup Date Range**: Filter by cage setup date
- **Export Filtered Results**: Download filtered cage data as CSV
- **Reset Filters**: One-click reset for all filter options

### ✅ 7. Gene Catalog Management
- **Gene Database**: Centralized gene and symbol management
- **Interactive Gene Table**: 
  - Checkbox selection for multiple genes
  - Click gene names to load into editor
  - Multi-select delete functionality
- **Quick Add Gene Flow**: 
  - Add genes on-the-fly during genotyping
  - Automatic continuation to genotype modal
  - Pre-selected genes in genotype forms
- **Symbol Management**: Comma-separated genotype symbols (e.g., `+/+, +/-, -/-, WT, flox/+`)

### ✅ 8. Enhanced Mouse Records Management
- **Comprehensive Filtering**:
  - Gene and genotype symbol filters
  - Date of birth range filtering
  - Sex, status, and cage filters
  - Free-text search across ID, notes, and genotype
  - Optional filters: Hide pups, Show only genotyped
- **Parent Tracking**: Mother_ID and Father_ID fields
- **Bulk Operations**: Multi-select and bulk delete functionality
- **Export Options**: Filtered data export in multiple formats

### ✅ 9. Improved Weaning and Pup Management
- **Strict Pup Definition**: Only Status = "Pup" counts as pups
- **Post-Wean Status Selection**: Choose Juvenile, Active, or Holding for weaned mice
- **Accurate Counts**: Proper pup vs total mouse counting
- **Flexible Weaning**: Support for multiple destination cages
- **Auto-Pick and Manual Selection**: Choose specific mice or auto-pick based on criteria

### ✅ 10. Breadcrumb Navigation
- **Hierarchical navigation**: Cage List > Cage Cx3-M1 > Mouse Cx3-M1-01
- **Clickable breadcrumbs** with pointer cursor
- **Visual feedback** on hover
- **Easy navigation** back to previous levels

### ✅ 11. Download & Export Features
- **Clear action labels**: "📥 Download Cage Data (Excel)"
- **Tooltips** explaining what each button does
- **Multiple export formats**: Excel, CSV, PDF, Print
- **Visual feedback** during export operations
- **Templated Exports**: mLIMS-compatible Excel exports

### ✅ 12. Advanced UI Features
- **Modal Management**: Smart modal flow for gene addition and genotyping
- **Selection Handling**: No more surprise modals when selecting rows
- **Click-to-Edit**: Click gene names to load into editor
- **Multi-Select Support**: Checkbox-based selection for bulk operations
- **Responsive Design**: Works on different screen sizes

### ✅ 13. Data Integrity & Validation
- **Accurate Pup Counting**: Strict definition prevents inflated counts
- **Proper Status Transitions**: Weaned mice get correct status
- **Filtered Data Operations**: Bulk operations work with filtered views
- **Duplicate Prevention**: Removed conflicting handlers and inputs

## 🏗️ System Architecture

### Core Modules
1. **Homepage** - Cage overview table with comprehensive filtering and keyboard shortcuts
2. **Cage View** - Detailed view of individual cages with mouse management
3. **Mouse Records** - Individual mouse tracking with advanced filtering
4. **Genes** - Gene catalog management with interactive table
5. **Breeding** - Breeding pair management with drag & drop
6. **Analytics** - Colony statistics and reporting
7. **Settings** - System configuration

### Data Structure
```r
# Cage Information (Enhanced Structure)
Cage_Name: "Breeding-05 (03/07)"
Genotype: "Cx3cr1(+/-)"
Sex: "M"
Num_Mice: 2
Num_Pups: 0
Has_Pups: "No"
Status: "Breeding"
Parent_Cage: ""
Genotyping_Done: TRUE
Genotype_Confirmed: TRUE
Setup_Date: "03/07"

# Mouse Information (Enhanced Structure)
Mouse_ID: "W1-01"
Cage: "Weaned-01"
Sex: "M"
DOB: "2024-03-15"
Genotype: "Cx3cr1(+/-)"
Status: "Juvenile"
Notes: "Weaned"
Mother_ID: "Breeding-05-01"
Father_ID: "Breeding-05-02"

# Gene Catalog Structure
Gene: "Cx3cr1"
Symbols: "+/+, +/-, -/-, WT, flox/+, flox/flox"
```

## 🚀 Installation & Setup

### Prerequisites
- R (version 4.0 or higher)
- RStudio (recommended)

### Required R Packages
```r
# Install required packages
install.packages(c(
  "shiny",
  "shinydashboard", 
  "DT",
  "dplyr",
  "readxl",
  "writexl",
  "shinyjs"
))
```

### Running the Application
1. **Open the project** in RStudio
2. **Install dependencies** if not already installed
3. **Run the application**:
   ```r
   # Method 1: From RStudio
   # Open app.R and click "Run App"
   
   # Method 2: From R console
   source("app.R")
   ```

## 🎮 Usage Guide

### Homepage Navigation
- **Click on cage names** in the table to view detailed information
- **Use keyboard shortcuts** for quick actions (Alt+N, Alt+D, Alt+P)
- **Press `/` key** to focus the search bar (when not typing)
- **Use comprehensive filters** to find specific cages
- **Export filtered results** for external analysis

### Cage Management
1. **View Cages**: Click on any cage name in the table to see details
2. **Add Cage**: Click "➕ New Cage" button or press `Alt + N`
3. **Search Cages**: Use the search bar or press `/` to focus
4. **Filter Cages**: Use the comprehensive filter panel
5. **Export Data**: Click "📥 Download All Data" or press `Alt + D`
6. **Print**: Click "🖨️ Print All" or press `Alt + P`

### Gene Catalog Management
1. **View Genes**: Navigate to the "🧬 Genes" tab
2. **Add Gene**: Use the gene form or quick-add during genotyping
3. **Edit Gene**: Click on gene names in the table to load into editor
4. **Delete Genes**: Select genes with checkboxes and click delete
5. **Manage Symbols**: Add comma-separated genotype symbols

### Mouse Records Management
1. **View Mice**: Navigate to the "🐭 Mouse Records" tab
2. **Filter Mice**: Use the comprehensive filter panel
3. **Add Mouse**: Click "➕ Add Mouse" button or press `Alt + M` (in cage view)
4. **Edit Mouse**: Press `Alt + Shift + E` to edit selected mouse
5. **Bulk Operations**: Select multiple mice for bulk delete
6. **Export Data**: Use table export buttons (Excel, CSV, PDF, Print)

### Genotyping Workflow
1. **Add Genotype**: Click "🧬 Add Genotype" in cage view
2. **No Genes?**: Automatically opens quick-add gene modal
3. **Select Gene**: Choose from gene catalog or add new one
4. **Choose Symbol**: Select from configured symbols or type custom
5. **Apply Scope**: Selected mice, all mice, or cage record only
6. **Save**: Genotype is applied and status updated

### Weaning Workflow
1. **Mark Weaned**: Click "📦 Wean / Split" in cage view
2. **Choose Status**: Select Juvenile, Active, or Holding for weaned mice
3. **Set Destinations**: Configure multiple destination cages
4. **Select Mice**: Choose specific mice or auto-pick
5. **Execute**: Mice are moved with correct status and counts updated

## 🎨 UI/UX Features

### Visual Design
- **Modern card-based layout**
- **Consistent color scheme** (blue primary, orange edit, red delete, green download)
- **Smooth animations** and transitions
- **Responsive design** for different screen sizes

### Accessibility Features
- **High contrast** color combinations
- **Clear typography** with good readability
- **Keyboard navigation** support
- **Screen reader** friendly structure
- **Focus indicators** for all interactive elements

### Interactive Elements
- **Hover effects** with color changes and shadows
- **Click feedback** with visual confirmation
- **Loading states** with wait cursors
- **Error states** with not-allowed cursors

## 📊 Data Management

### Import/Export Capabilities
- **Excel files** (.xlsx, .xls) with mLIMS templates
- **CSV files** (.csv)
- **PDF reports** for printing
- **Direct printing** functionality
- **Filtered data exports**

### Data Validation
- **Input validation** for all forms
- **Error handling** for invalid data
- **Confirmation dialogs** for destructive actions
- **Auto-save** functionality

## 🔧 Customization

### Styling
The application uses CSS classes that can be easily customized:
```css
.clickable { cursor: pointer; }
.cage-card { /* Cage styling */ }
.action-btn { /* Button styling */ }
.breadcrumb { /* Navigation styling */ }
```

### Adding New Features
1. **New tabs**: Add to `sidebarMenu` in UI
2. **New functionality**: Add to server logic
3. **New data types**: Extend data structures
4. **New exports**: Add to export functions

## 🐛 Troubleshooting

### Common Issues
1. **Package not found**: Install missing packages
2. **Port conflicts**: Change port in `runApp()`
3. **Data not loading**: Check file paths and permissions
4. **UI not responsive**: Check browser compatibility

### Performance Tips
- **Limit data size** for large colonies
- **Use pagination** for large tables
- **Optimize queries** for database operations
- **Cache frequently used data**

## 📝 Development Notes

### Code Structure
- **Modular design** for easy maintenance
- **Reactive programming** for dynamic updates
- **Event-driven architecture** for user interactions
- **Separation of concerns** between UI and logic

### Best Practices
- **Consistent naming** conventions
- **Error handling** throughout
- **User feedback** for all actions
- **Accessibility** considerations
- **Performance** optimization

## 🤝 Contributing

### Development Workflow
1. **Fork the repository**
2. **Create feature branch**
3. **Make changes** with proper testing
4. **Submit pull request**
5. **Code review** and approval

### Testing Guidelines
- **Test all interactive elements**
- **Verify accessibility features**
- **Check responsive design**
- **Validate data operations**

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🆘 Support

For support and questions:
- **Documentation**: Check this README
- **Issues**: Report bugs via GitHub issues
- **Features**: Request new features via GitHub issues
- **Email**: Contact the development team

---

**Version**: 2.0.0  
**Last Updated**: January 2024  
**Compatibility**: R 4.0+, Shiny 1.7+
