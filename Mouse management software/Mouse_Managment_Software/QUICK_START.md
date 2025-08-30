# Quick Start Guide - Zijie's Great Mouse Management Software

## 🚀 Get Started in 5 Minutes

### Step 1: Install R and RStudio
1. **Download R** from [r-project.org](https://www.r-project.org/)
2. **Download RStudio** from [rstudio.com](https://www.rstudio.com/)
3. **Install both** following the installation guides

### Step 2: Install Required Packages
```r
# Run this in R console or RStudio
source("install_packages.R")
```

### Step 3: Run the Application
```r
# Method 1: From RStudio
# Open app.R and click "Run App" button

# Method 2: From R console
source("app.R")
```

### Step 4: Explore the Features
- **Homepage**: Cage overview table with comprehensive filtering and shortcuts
- **Cage View**: Detailed view of individual cages with mouse management
- **Mouse Records**: Advanced filtering and bulk operations
- **Genes**: Gene catalog management with interactive table
- **Breeding**: Manage breeding pairs with drag & drop
- **Analytics**: Generate reports and statistics
- **Settings**: Configure the system

## 🎯 Key Features to Try

### 1. Smart Keyboard Shortcuts
- **Type freely**: Use `/` in input fields without interference (e.g., type `flox/+`)
- **Quick actions**: Alt+N (New Cage), Alt+D (Download), Alt+M (Add Mouse)
- **Search focus**: Press `/` when not typing to focus search bar
- **Navigation**: Esc or Alt+H to go back to home

### 2. Comprehensive Filtering System
- **Cage Filtering**: Filter by gene, genotype, sex, status, has pups, genotyping status, setup date
- **Mouse Filtering**: Filter by gene, genotype, DOB range, sex, status, cage, search text
- **Export filtered results**: Download only the data you need
- **Reset filters**: One-click reset for all filter options

### 3. Gene Catalog Management
- **Add genes**: Use the gene form or quick-add during genotyping
- **Interactive table**: Click gene names to load into editor, use checkboxes for multi-select
- **Symbol management**: Add comma-separated genotype symbols (e.g., `+/+, +/-, -/-, WT, flox/+`)
- **Smart workflow**: No genes? Automatically opens quick-add modal

### 4. Enhanced Mouse Records
- **Advanced filtering**: Find mice by any combination of criteria
- **Bulk operations**: Select multiple mice for bulk delete
- **Parent tracking**: Mother_ID and Father_ID fields
- **Optional filters**: Hide pups by default, show only genotyped

### 5. Improved Weaning Workflow
- **Post-wean status**: Choose Juvenile, Active, or Holding for weaned mice
- **Accurate counts**: Proper pup vs total mouse counting
- **Multiple destinations**: Support for splitting into multiple cages
- **Flexible selection**: Choose specific mice or auto-pick

### 6. Click Targets & Hover Effects
- **Hover over cage names** in the table to see them highlight
- **Click on cage names** to view detailed cage information
- **Notice cursor changes** to pointer (🖱️) on clickable elements
- **No surprise modals**: Select rows without opening detail modals

### 7. Navigation
- **Use Tab key** to navigate through elements
- **Click breadcrumbs** to navigate back: Home > Cage Cx3-M1 > Mouse M1-01
- **Use Enter/Space** to activate buttons
- **Press `/` key** to focus search bar (when not typing)
- **Use keyboard shortcuts** for quick actions

### 8. Export Features
- **Click download buttons** (📥) on homepage for all data
- **Use keyboard shortcuts**: Alt+D (Download All), Alt+E (Export Cage)
- **Export filtered results**: Download only the data you need
- **Use table export** in Mouse Records tab
- **Generate reports** in Analytics tab

### 9. Interactive Elements
- **Drag & drop area** in Breeding tab (cursor changes to move ↕️)
- **Tooltips** appear on hover over icons
- **Loading states** show wait cursor (⏳) during operations
- **Parent cage links** clickable in cage table
- **Quick action buttons** in cage view for common tasks
- **Gene table interactions**: Click genes to edit, use checkboxes for bulk operations

## 🎨 UI/UX Features Demonstrated

### Cursor Behavior Guide
| Element | Cursor | Action |
|---------|--------|--------|
| Cage names | Pointer (🖱️) | Click to view details |
| Mouse IDs | Pointer (🖱️) | Click to view mouse info |
| Gene names | Pointer (🖱️) | Click to load into editor |
| Download buttons | Pointer (🖱️) | Click to export data |
| Drag areas | Move (↕️) | Drag to create pairs |
| Loading | Wait (⏳) | Operation in progress |
| Disabled | Not allowed (🚫) | Action not available |

### Visual Feedback
- **Hover highlights**: Elements change color and lift
- **Click feedback**: Visual confirmation of actions
- **Smooth transitions**: All animations are fluid
- **Color coding**: Blue (view), Orange (edit), Red (delete), Green (download)

## 📊 Sample Data Included

The application comes with sample data including:
- **6 cages** including breeding and weaned cages with parent relationships
- **8 mice** including weaned pups with unknown genotypes
- **2 genes** (Cx3cr1, Tmem119) with genotype symbols
- **Breeding status** tracking with multiple status types
- **Genotyping status** with pending and confirmed states
- **Parent-child relationships** between breeding and weaned cages
- **Enhanced data structure** with Mother_ID, Father_ID, and detailed fields

## 🔧 Customization

### Adding Your Data
1. **Replace sample data** in `sample_data.R`
2. **Modify data structures** to match your needs
3. **Update UI elements** in `app.R` if needed
4. **Add your genes** to the gene catalog

### Styling Changes
- **Edit CSS** in the `tags$style` section of `app.R`
- **Modify colors** by changing hex values
- **Adjust animations** by modifying transition properties

## 🐛 Troubleshooting

### Common Issues
1. **"Package not found"**
   - Run `source("install_packages.R")` again
   - Check internet connection

2. **"Port already in use"**
   - Close other R applications
   - Change port in `runApp(port = 8080)`

3. **"UI not loading"**
   - Check browser compatibility
   - Try refreshing the page

4. **"Keyboard shortcuts not working"**
   - Make sure you're not typing in an input field
   - Check that the app has focus

### Getting Help
- **Check README.md** for detailed documentation
- **Review sample_data.R** for data structure examples
- **Examine app.R** for code examples

## 📈 Next Steps

### Advanced Features to Add
1. **Database integration** (SQLite, PostgreSQL)
2. **User authentication** and permissions
3. **Advanced analytics** and visualizations
4. **Email notifications** for health checks
5. **Mobile-responsive design**

### Data Management
1. **Backup systems** for data safety
2. **Import/export** from existing systems
3. **Data validation** and error checking
4. **Audit trails** for changes

## 🎉 You're Ready!

Your mouse colony management software is now running with all the advanced features:

✅ **Smart keyboard shortcuts** that don't interfere with typing  
✅ **Comprehensive filtering** for cages and mice  
✅ **Gene catalog management** with interactive table  
✅ **Enhanced mouse records** with advanced filtering  
✅ **Improved weaning workflow** with accurate counts  
✅ **Click-to-edit functionality** for genes and mice  
✅ **Bulk operations** for efficient data management  
✅ **Modal management** with smart workflow  
✅ **Data integrity** with proper validation  
✅ **Export capabilities** for filtered data  
✅ **Responsive design** with modern UI/UX  

Start managing your mouse colony efficiently and intuitively with the most advanced features available!

## 🚀 Quick Workflow Examples

### Adding a New Gene and Genotyping
1. Click "🧬 Add Genotype" in any cage
2. If no genes exist, the quick-add gene modal opens
3. Enter gene name (e.g., "MyGene") and symbols (e.g., "+/+, +/-, -/-, WT")
4. Click "Save & Continue" - genotype modal opens with gene pre-selected
5. Choose genotype symbol and apply to mice or cage
6. Done!

### Filtering and Exporting Data
1. Use the filter panel on Home or Cage View
2. Select gene, status, has pups, etc.
3. Click "📥 Export filtered" to download only filtered data
4. Use "Reset filters" to clear all filters

### Managing Mouse Records
1. Go to "🐭 Mouse Records" tab
2. Use filters to find specific mice
3. Select multiple mice with checkboxes
4. Use bulk delete or export filtered data
5. Click mouse IDs to view detailed information

### Weaning Mice
1. Click "📦 Wean / Split" in cage view
2. Choose post-wean status (Juvenile, Active, or Holding)
3. Set up destination cages and mouse counts
4. Select specific mice or let system auto-pick
5. Execute - mice are moved with correct status and counts updated
