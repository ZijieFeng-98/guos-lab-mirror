# Mouse Colony Management Software (mLIMS)
# Main Shiny Application with Keyboard Shortcuts

# Allow up to ~200 MB uploads (tweak as you like)
options(shiny.maxRequestSize = 200 * 1024^2)

library(shiny)
library(shinydashboard)
library(DT)
library(dplyr)
library(readxl)
library(writexl)
library(shinyjs)

# Database and backup packages (optional - will work without DB)
# library(DBI)
# library(RPostgres)
# library(pool)
# library(pins)        # for versioned backups to S3/GCS/etc.
# library(jsonlite)    # tiny helper for safe text fields if needed

# UI Definition
ui <- dashboardPage(
  dashboardHeader(title = "Zijie's Great Mouse Management Software", titleWidth = 350),
  
  dashboardSidebar(
    width = 300,
          sidebarMenu(
        id = "sidebar",
        menuItem("🏠 Home", tabName = "dashboard", icon = icon("home")),
        menuItem("🔎 Cage View", tabName = "cage_view", icon = icon("binoculars")),
        menuItem("🐭 Mouse Records", tabName = "mice", icon = icon("mouse")),
        menuItem("📊 Breeding", tabName = "breeding", icon = icon("heart")),
        menuItem("📈 Analytics", tabName = "analytics", icon = icon("chart-line")),
        menuItem("🧬 Genes", tabName = "genes", icon = icon("dna")),
        menuItem("⚙️ Settings", tabName = "settings", icon = icon("cog"))
      )
  ),
  
  dashboardBody(
    useShinyjs(),
    # Hidden heartbeat (keeps session lightly ticking)
    div(style = "display:none;", textOutput("heartbeat")),
    # Offline banner
    div(id = "offlineBanner", "Connection lost… trying to reconnect automatically."),
    tags$head(
      tags$title("Zijie's Great Mouse Management Software"),
      tags$style(HTML("
        #offlineBanner {
          display:none; position:fixed; bottom:16px; left:16px; right:16px;
          background:#ffebee; color:#b71c1c; border:1px solid #ffcdd2;
          padding:10px 14px; border-radius:8px; z-index:9999; font-weight:600;
        }
      ")),
      tags$script(HTML("
        // --- Resume state (tab + current cage) across reloads ---
        Shiny.addCustomMessageHandler('rememberState', function(x){
          try { localStorage.setItem('mlims_state', JSON.stringify(x || {})); } catch(e){}
        });

        function restoreState(){
          try {
            var raw = localStorage.getItem('mlims_state'); if (!raw) return;
            var x = JSON.parse(raw);
            if (x && x.tab)  Shiny.setInputValue('restore_tab',  x.tab,  {priority:'event'});
            if (x && x.cage) Shiny.setInputValue('restore_cage', x.cage, {priority:'event'});
          } catch(e){}
        }

        // Restore after page load or when connection comes back
        document.addEventListener('DOMContentLoaded', restoreState);
        $(document).on('shiny:connected', restoreState);

        // --- Auto reload after disconnect (fallback if Shiny can't reconnect) ---
        var autoReloadTimer = null;
        $(document).on('shiny:disconnected', function(){
          // Try a gentle auto-reload after a few seconds.
          // (User can still click the default 'Reload' button.)
          if (autoReloadTimer) clearTimeout(autoReloadTimer);
          autoReloadTimer = setTimeout(function(){ location.reload(); }, 4000);
        });

        // Show offline banner
        $(document).on('shiny:disconnected', function(){ $('#offlineBanner').fadeIn(150); });
        $(document).on('shiny:connected',    function(){ $('#offlineBanner').fadeOut(150); });

        // Keyboard shortcuts for mLIMS
        document.addEventListener('keydown', function(event) {
          const ae = document.activeElement;
          const tag = ae && ae.tagName ? ae.tagName.toLowerCase() : '';
          const isEditable = tag === 'input' || tag === 'textarea' || tag === 'select' || (ae && ae.isContentEditable);

          // Allow Esc to work anywhere; ignore other shortcuts while typing in fields
          if (isEditable && event.key !== 'Escape') return;

          // Alt + N: New Cage
          if (event.altKey && event.key === 'n') {
            event.preventDefault();
            Shiny.setInputValue('shortcut_new_cage', Math.random(), {priority: 'event'});
          }

          // Alt + D: Download All Data
          if (event.altKey && event.key === 'd') {
            event.preventDefault();
            Shiny.setInputValue('shortcut_download_all', Math.random(), {priority: 'event'});
          }

          // Alt + P: Print All
          if (event.altKey && event.key === 'p' && !event.shiftKey) {
            event.preventDefault();
            Shiny.setInputValue('shortcut_print_all', Math.random(), {priority: 'event'});
          }

          // Alt + Shift + P: Print Current Cage
          if (event.altKey && event.shiftKey && event.key === 'p') {
            event.preventDefault();
            Shiny.setInputValue('shortcut_print_cage', Math.random(), {priority: 'event'});
          }

          // Alt + M: Add Mouse
          if (event.altKey && event.key === 'm') {
            event.preventDefault();
            Shiny.setInputValue('shortcut_add_mouse', Math.random(), {priority: 'event'});
          }

          // Alt + E: Export Current Cage
          if (event.altKey && event.key === 'e') {
            event.preventDefault();
            Shiny.setInputValue('shortcut_export_cage', Math.random(), {priority: 'event'});
          }

          // Alt + Shift + E: Edit Mouse Record
          if (event.altKey && event.shiftKey && event.key === 'e') {
            event.preventDefault();
            Shiny.setInputValue('shortcut_edit_mouse', Math.random(), {priority: 'event'});
          }

          // Alt + U: Upload Data
          if (event.altKey && event.key === 'u') {
            event.preventDefault();
            Shiny.setInputValue('shortcut_upload', Math.random(), {priority: 'event'});
          }

          // Alt + H or Esc: Back to Home or Cage List
          if ((event.altKey && event.key === 'h') || event.key === 'Escape') {
            event.preventDefault();
            Shiny.setInputValue('shortcut_back_home', Math.random(), {priority: 'event'});
          }

          // '/' focuses search ONLY when not editing a field
          if (event.key === '/' && !event.ctrlKey && !event.altKey && !isEditable) {
            event.preventDefault();
            var cand = ['#cage_search', '#cage_search_cv'];
            for (var i=0;i<cand.length;i++){
              var el = document.querySelector(cand[i]);
              if (el && el.offsetParent !== null) { el.focus(); break; } // visible one
            }
          }
        });
      ")),
      tags$style(HTML("
        /* Global cursor and hover styles */
        .clickable {
          cursor: pointer !important;
          transition: all 0.2s ease;
        }
        
        .clickable:hover {
          background-color: #e3f2fd !important;
          border-color: #2196f3 !important;
          transform: translateY(-1px);
          box-shadow: 0 2px 8px rgba(33, 150, 243, 0.3);
        }
        
        /* Keyboard shortcut hints */
        .shortcut-hint {
          font-size: 11px;
          color: #666;
          margin-left: 5px;
        }
        
        /* Search bar styling */
        .search-container {
          margin-bottom: 20px;
        }
        
        .search-container .form-control {
          width: 100%;
          padding: 10px;
          border: 2px solid #ddd;
          border-radius: 4px;
          font-size: 14px;
        }
        
        .search-container .form-control:focus {
          border-color: #2196f3;
          outline: none;
        }
        
        .cage-card {
          border: 2px solid #e0e0e0;
          border-radius: 8px;
          padding: 15px;
          margin: 10px 0;
          background: white;
          cursor: pointer;
          transition: all 0.3s ease;
        }
        
        .cage-card:hover {
          border-color: #2196f3;
          background-color: #f8f9ff;
          transform: translateY(-2px);
          box-shadow: 0 4px 12px rgba(33, 150, 243, 0.2);
        }
        
        .action-btn {
          cursor: pointer;
          padding: 8px 16px;
          border-radius: 4px;
          border: none;
          font-weight: 500;
          transition: all 0.2s ease;
        }
        
        .btn-view {
          background-color: #2196f3;
          color: white;
        }
        
        .btn-view:hover {
          background-color: #1976d2;
          transform: translateY(-1px);
        }
        
        .btn-edit {
          background-color: #ff9800;
          color: white;
        }
        
        .btn-edit:hover {
          background-color: #f57c00;
          transform: translateY(-1px);
        }
        
        .btn-delete {
          background-color: #f44336;
          color: white;
        }
        
        .btn-delete:hover {
          background-color: #d32f2f;
          transform: translateY(-1px);
        }
        
        .btn-download {
          background-color: #4caf50;
          color: white;
        }
        
        .btn-download:hover {
          background-color: #388e3c;
          transform: translateY(-1px);
        }
        
        .btn-danger {
          background-color: #f44336;
          color: white;
        }
        
        .btn-danger:hover {
          background-color: #d32f2f;
          transform: translateY(-1px);
        }
        
        .breadcrumb {
          background-color: #f5f5f5;
          padding: 10px 15px;
          border-radius: 4px;
          margin-bottom: 20px;
        }
        
        .breadcrumb-item {
          color: #2196f3;
          cursor: pointer;
          text-decoration: none;
        }
        
        .breadcrumb-item:hover {
          text-decoration: underline;
          color: #1976d2;
        }
        
        .breadcrumb-separator {
          color: #666;
          margin: 0 8px;
        }
        
        .mouse-id {
          cursor: pointer;
          color: #2196f3;
          font-weight: 500;
        }
        
        .mouse-id:hover {
          text-decoration: underline;
          color: #1976d2;
        }
        
        .tooltip {
          position: relative;
          display: inline-block;
        }
        
        .tooltip .tooltiptext {
          visibility: hidden;
          width: 200px;
          background-color: #333;
          color: white;
          text-align: center;
          border-radius: 6px;
          padding: 5px;
          position: absolute;
          z-index: 1;
          bottom: 125%;
          left: 50%;
          margin-left: -100px;
          opacity: 0;
          transition: opacity 0.3s;
        }
        
        .tooltip:hover .tooltiptext {
          visibility: visible;
          opacity: 1;
        }
        
        .loading {
          cursor: wait !important;
        }
        
        .disabled {
          cursor: not-allowed !important;
          opacity: 0.6;
        }
        
        .drag-area {
          cursor: move;
          border: 2px dashed #ccc;
          padding: 20px;
          text-align: center;
          border-radius: 8px;
        }
        
        .drag-area:hover {
          border-color: #2196f3;
          background-color: #f8f9ff;
        }
        
        /* Right-align header title */
        .main-header .logo { text-align: right; }
      "))
    ),
    
    tabItems(
      # Homepage: Cage Overview Tab
      tabItem(tabName = "dashboard",
        fluidRow(
          column(12,
            div(class = "breadcrumb",
              span("🏠 Home", class = "breadcrumb-item",
                 onclick = "Shiny.setInputValue('breadcrumb_click', 'home', {priority:'event'})")
            )
          )
        ),
        fluidRow(
          column(12,
            div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;",
              h2("Cage Overview"),
              div(
                actionButton("new_cage_home", "➕ New Cage", class = "action-btn btn-view",
                           onclick = "Shiny.setInputValue('shortcut_new_cage', 1, {priority: 'event'})"),
                actionButton("open_upload", "⬆️ Upload Data", class = "action-btn btn-view",
                           onclick = "Shiny.setInputValue('shortcut_upload', 1, {priority: 'event'})"),
                downloadButton("dl_all", "📥 Download All Data", class = "action-btn btn-download"),
                downloadButton("dl_cages_mlims_xlsx", "📥 Cages (Excel)", class = "action-btn btn-download"),
                downloadButton("dl_labmate_all", "📥 Mice (Labmate format)", class = "action-btn btn-download"),
                actionButton("print_all_home", "🖨️ Print All", class = "action-btn btn-download",
                           onclick = "Shiny.setInputValue('shortcut_print_all', 1, {priority: 'event'})"),
                actionButton("delete_selected_cages_dashboard", "🗑️ Delete Selected", class = "action-btn btn-danger")
              )
            )
          )
        ),
        fluidRow(
          column(12,
            div(class = "search-container",
              textInput("cage_search", "Search Cages (Press '/' to focus)", 
                       placeholder = "Type to search cages...")
            )
          )
        ),
        fluidRow(
          column(12,
            div(style="margin:8px 0 16px 0; font-size:12px; color:#666;",
                HTML("Tip: Press <b>Alt+U</b> to open the Import dialog."))
          )
        ),
        fluidRow(
          column(12,
            div(class="box", style="padding:12px; border:1px solid #eee; border-radius:8px; background:#fafafa;",
              h4("Filters"),
              fluidRow(
                column(3, uiOutput("cage_filter_gene_ui")),
                column(3, uiOutput("cage_filter_symbol_ui")),
                column(3, selectInput("cage_filter_sex", "Sex", c("Any","M","F","Mixed"), "Any")),
                column(3, selectInput("cage_filter_status", "Status",
                                      c("Any","Breeding","Holding","Experiment","Active","Retired"), "Any"))
              ),
              fluidRow(
                column(3, selectInput("cage_filter_has_pups", "Has pups", c("Any","Yes","No"), "Any")),
                column(3, selectInput("cage_filter_genotyped", "Genotyping done", c("Any","Yes","No"), "Any")),
                column(3, selectInput("cage_filter_confirmed", "Genotype confirmed", c("Any","Yes","No"), "Any")),
                column(3, dateRangeInput("cage_filter_setup", "Setup date", start = Sys.Date()-60, end = Sys.Date()))
              ),
              div(style="display:flex; gap:8px; margin-top:8px;",
                actionButton("cage_filters_reset", "Reset filters", class="action-btn btn-edit"),
                downloadButton("dl_filtered_cages", "📥 Export filtered", class="action-btn btn-download"),
                div(style="flex:1"),
                checkboxInput("filter_genotyping", "Show only cages with pending genotyping", value = FALSE) # keep your quick toggle
              )
            )
          )
        ),
        fluidRow(
          column(12,
            div(style = "margin-bottom: 20px;",
              span("💡 Keyboard shortcuts: Alt+N (New Cage), Alt+D (Download), Alt+P (Print), / (Search)", 
                  class = "shortcut-hint")
            )
          )
        ),
        fluidRow(
          column(12,
            DTOutput("cage_overview_table")
          )
        )
      ),
      
      # Cage View Tab (Dynamic - shows list or detail view)
      tabItem(tabName = "cage_view",
        uiOutput("cage_view_ui")
      ),
      
      # Mouse Records Tab
      tabItem(tabName = "mice",
        fluidRow(
          column(12,
            div(class = "breadcrumb",
              span("🐭 Mouse Records", class = "breadcrumb-item")
            )
          )
        ),
        fluidRow(
          column(12,
            div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;",
              h2("Mouse Records"),
              div(
                actionButton("add_mouse", "➕ Add Mouse", class = "action-btn btn-view"),
                actionButton("import_mice", "📥 Import Mice", class = "action-btn btn-download"),
                downloadButton("dl_mice_all_mlims_xlsx", "📤 All Mice (mLIMS Excel)", class = "action-btn btn-download"),
                downloadButton("dl_mice_filtered_mlims_xlsx", "📤 Filtered (mLIMS Excel)", class = "action-btn btn-download"),
                downloadButton("dl_labmate_xlsx", "📤 Mice (Labmate format)", class = "action-btn btn-download"),
                actionButton("delete_selected_mice_main", "🗑️ Delete Selected", class = "action-btn btn-danger")
              )
            )
          )
        ),
        fluidRow(
          column(12,
            div(class="box", style="padding:12px; border:1px solid #eee; border-radius:8px; background:#fafafa;",
              h4("Filters"),
              fluidRow(
                column(3, uiOutput("mouse_filter_gene_ui")),
                column(3, uiOutput("mouse_filter_symbol_ui")),
                column(3, dateRangeInput("mouse_filter_dob", "DOB range", start = Sys.Date() - 90, end = Sys.Date())),
                column(3, selectInput("mouse_filter_sex", "Sex", c("Any","M","F"), "Any"))
              ),
              fluidRow(
                column(3, selectInput("mouse_filter_status", "Status",
                                      c("Any","Pup","Juvenile","Active","Holding","Experiment","Retired"), "Any")),
                column(3, selectizeInput("mouse_filter_cage", "Cage", choices = NULL, multiple = TRUE,
                                         options = list(placeholder="Any cage"))),
                column(6, textInput("mouse_search_text", "Search (ID / notes / genotype)", ""))
              ),
              div(style="display:flex; gap:8px; margin-top:8px;",
                actionButton("mouse_filters_reset", "Reset filters", class="action-btn btn-edit"),
                downloadButton("dl_filtered_mice", "📥 Export filtered", class="action-btn btn-download")
              ),
              div(style="display:flex; gap:16px; margin-top:8px; font-size:12px;",
                checkboxInput("hide_pups_mice", "Hide pups (pre-wean)", TRUE),
                checkboxInput("only_genotyped_mice", "Show only genotyped", FALSE)
              )
            )
          )
        ),
        fluidRow(
          column(12,
            DTOutput("mice_table")
          )
        )
      ),
      
      # Breeding Tab
      tabItem(tabName = "breeding",
        fluidRow(
          column(12,
            div(class = "breadcrumb",
              span("📊 Breeding", class = "breadcrumb-item")
            )
          )
        ),
        fluidRow(
          column(12,
            h2("Breeding Management"),
            div(class = "drag-area",
              h4("🔄 Drag & Drop Breeding Pairs"),
              p("Drag mice here to create breeding pairs"),
              p("Cursor changes to move (↕️) when dragging")
            )
          )
        )
      ),
      
      # Analytics Tab
      tabItem(tabName = "analytics",
        fluidRow(
          column(12,
            div(class = "breadcrumb",
              span("📈 Analytics", class = "breadcrumb-item")
            )
          )
        ),
        fluidRow(
          column(12,
            h2("Colony Analytics"),
            div(style = "margin-bottom: 16px;",
              actionButton("generate_report", "📊 Generate Report", class = "action-btn btn-view"),
              actionButton("print_view", "🖨️ Print View", class = "action-btn btn-download")
            )
          )
        ),

        # --- Daily Log panel ---
        fluidRow(
          column(12,
            div(class="box", style="padding:12px; border:1px solid #eee; border-radius:8px; background:#fafafa;",
              h3("📝 Daily Log"),
              fluidRow(
                column(4, dateInput("log_date", "Date", value = Sys.Date())),
                column(4, uiOutput("log_action_filter_ui")),
                column(4, div(style="margin-top:25px;",
                              downloadButton("dl_daily_log_html", "⬇️ Download Printable HTML", class="action-btn btn-download"),
                              downloadButton("dl_audit_csv", "⬇️ Download Audit CSV", class="action-btn btn-download")
                      ))
              ),
              div(style="color:#666; font-size:12px; margin-top:4px;",
                  "Tip: open the HTML and press Ctrl/Cmd+P to print or save as PDF."),
              hr(),
              uiOutput("daily_log_summary_ui"),
              DTOutput("daily_log_table")
            )
          )
        )
      ),
      
      # Genes Tab
      tabItem(tabName = "genes",
        fluidRow(
          column(12,
            div(class = "breadcrumb",
              span("🧬 Genes", class = "breadcrumb-item")
            )
          )
        ),
        fluidRow(
          column(12,
            h2("Gene Catalog"),
            div(style="margin:10px 0 16px 0; color:#666; font-size:13px;",
                "Add a gene and the allowed genotype symbols (comma-separated). ",
                "These symbols will be the choices when you assign genotypes."
            ),
            fluidRow(
              column(4, textInput("gene_name", "Gene name:", placeholder = "e.g., Cx3cr1")),
              column(6, textInput("gene_symbols", "Allowed symbols (comma-separated):",
                                  placeholder = "e.g., +/+, +/-, -/-, WT, flox/+, flox/flox")),
              column(2, div(style="margin-top:24px;",
                            actionButton("add_gene_btn", "➕ Add / Update Gene", class="action-btn btn-view")))
            ),
            div(style="display:flex; gap:8px; align-items:center; margin:6px 0 12px 0;",
              actionButton("delete_gene_btn", "🗑️ Delete Selected", class="action-btn btn-danger")
            ),
            DTOutput("gene_table")
          )
        )
      ),
      
      # Settings Tab
      tabItem(tabName = "settings",
        fluidRow(
          column(12,
            div(class = "breadcrumb",
              span("⚙️ Settings", class = "breadcrumb-item")
            )
          )
        ),
        fluidRow(
          column(12,
            h2("System Settings"),
            div(style = "margin-bottom: 20px;",
              actionButton("save_settings", "💾 Save Settings", class = "action-btn btn-view"),
              actionButton("reset_settings", "🔄 Reset to Default", class = "action-btn btn-edit")
            )
          )
        ),
        
        # --- Backup & Restore Section ---
        fluidRow(
          column(12,
            div(class="box", style="padding:12px; border:1px solid #eee; border-radius:8px; background:#fafafa;",
              h3("☁️ Backup & Restore"),
              p("Create versioned backups of your data to S3 cloud storage."),
              div(style = "margin-bottom: 16px;",
                actionButton("backup_now", "☁️ Backup Now", class="action-btn btn-download"),
                actionButton("restore_latest", "🔁 Restore Latest", class="action-btn btn-edit")
              ),
              div(style="color:#666; font-size:12px;",
                  "Backups include all cages, mice, genes, and audit logs. Each backup is versioned and immutable.")
            )
          )
        ),
        
        # --- Database Status ---
        fluidRow(
          column(12,
            div(class="box", style="padding:12px; border:1px solid #eee; border-radius:8px; background:#fafafa;",
              h3("🗄️ Database Status"),
              uiOutput("db_status_ui")
            )
          )
        )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # Reactive values for app state
  current_cage <- reactiveVal(NULL)
  current_mouse <- reactiveVal(NULL)
  
  # Initial data frames (defined first, before reactiveVal)
  initial_cages <- data.frame(
    Cage_Name = c("Cx3-M1", "Tmem-M1", "Cx3-M2", "Tmem-M2"),
    Genotype = c("Cx3cr1", "Tmem119", "Cx3cr1", "Tmem119"),
    Sex = c("M", "M", "F", "F"),
    Num_Mice = c(1, 1, 2, 1),
    Num_Pups = c(0, 0, 0, 0),
    Status = c("Holding", "Experiment", "Active", "Active"),
    Parent_Cage = c("", "", "", ""),
    Genotyping_Done = c(TRUE, TRUE, FALSE, TRUE),
    Genotype_Confirmed = c(TRUE, TRUE, FALSE, FALSE),
    Setup_Date = c("", "", "", ""),
    Notes = c("", "", "", ""),
    stringsAsFactors = FALSE
  )
  
  initial_mice <- data.frame(
    Mouse_ID = c("M1-01", "M1-02", "M2-01", "M2-02"),
    Cage = c("Cx3-M1", "Cx3-M1", "Tmem-M1", "Tmem-M1"),
    Sex = c("M", "F", "M", "F"),
    DOB = c("7/25", "7/25", "8/1", "8/1"),
    Genotype = c("+/+", "+/+", "+/+", "+/+"),
    Status = c("Active", "Active", "Active", "Active"),
    Notes = c("Transferred", "Transferred", "New", "New"),
    Mother_ID = NA_character_,
    Father_ID = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Reactive data storage initialized with initial data
  cage_data_reactive <- reactiveVal(initial_cages)
  mice_data_reactive <- reactiveVal(initial_mice)
  
  # Load from database if available, otherwise use initial data
  is_db_empty <- function(df) is.null(df) || !nrow(df)
  
  # Initialize database connection and load data
  # (moved after all reactiveVals are defined)
  
  # Always keep counts fresh (also runs once on app start)
  observeEvent(mice_data_reactive(), {
    recalc_cage_counts()
  }, ignoreInit = FALSE)
  
  # Holds the template header (column names) if present
  mlims_cols <- reactiveVal(NULL)
  
  # Define a cage template (minimal, app-friendly)
  cage_template_cols <- c(
    "Cage_Name","Genotype","Sex","Num_Mice","Num_Pups","Status",
    "Parent_Cage","Setup_Date","Genotyping_Done","Genotype_Confirmed","Notes"
  )
  
  # Upload data reactive values
  uploaded_df <- reactiveVal(NULL)
  upload_issues <- reactiveVal(list(errors = character(), warnings = character()))
  
  # Save last state for quick undo
  prev_cages <- reactiveVal(NULL)
  prev_mice  <- reactiveVal(NULL)
  
  # Gene catalog: Gene + Symbols (comma-separated)
  gene_catalog <- reactiveVal(
    data.frame(
      Gene    = c("Cx3cr1","Tmem119"),
      Symbols = c("+/+,+/-,-/-,WT", "+/+,+/-,-/-,WT"),
      stringsAsFactors = FALSE
    )
  )
  last_gene <- reactiveVal(NULL)  # remember last used gene in genotype dialog
  
  # --- Daily Log / Audit Trail -----------------------------------------------
  audit_log <- reactiveVal(
    data.frame(
      Time      = as.POSIXct(character()),
      Action    = character(),
      Cage      = character(),
      Mouse_ID  = character(),
      Count     = integer(),
      Details   = character(),
      stringsAsFactors = FALSE
    )
  )

  log_event <- function(action, cage = NA_character_, mouse = NA_character_,
                        details = "", count = NA_integer_) {
    df <- audit_log()
    row <- data.frame(
      Time     = Sys.time(),
      Action   = as.character(action %||% ""),
      Cage     = as.character(cage %||% NA_character_),
      Mouse_ID = as.character(mouse %||% NA_character_),
      Count    = as.integer(count %||% NA_integer_),
      Details  = as.character(details %||% ""),
      stringsAsFactors = FALSE
    )
    audit_log(rbind(df, row))
  }
  
  # ---------- Persistent Storage (Postgres) ----------
  make_db_pool <- function() {
    pool::dbPool(
      drv      = RPostgres::Postgres(),
      host     = Sys.getenv("DB_HOST"),
      port     = as.integer(Sys.getenv("DB_PORT", "5432")),
      dbname   = Sys.getenv("DB_NAME"),
      user     = Sys.getenv("DB_USER"),
      password = Sys.getenv("DB_PASSWORD"),
      sslmode  = Sys.getenv("DB_SSLMODE", "require")
    )
  }

  # Try to create pool, fall back gracefully if DB not configured
  pool <- tryCatch({
    if (nzchar(Sys.getenv("DB_HOST"))) {
      make_db_pool()
    } else {
      NULL
    }
  }, error = function(e) {
    message("Database not configured or unavailable: ", e$message)
    NULL
  })

  init_db <- function(conn) {
    if (is.null(conn)) return()
    
    DBI::dbExecute(conn, "
      CREATE TABLE IF NOT EXISTS cages (
        cage_name           TEXT PRIMARY KEY,
        genotype            TEXT,
        sex                 TEXT,
        num_mice            INT,
        num_pups            INT,
        status              TEXT,
        parent_cage         TEXT,
        genotyping_done     BOOLEAN,
        genotype_confirmed  BOOLEAN,
        setup_date          TEXT,
        notes               TEXT,
        updated_at          TIMESTAMPTZ DEFAULT now()
      );")

    DBI::dbExecute(conn, "
      CREATE TABLE IF NOT EXISTS mice (
        mouse_id         TEXT PRIMARY KEY,
        cage             TEXT REFERENCES cages(cage_name) ON DELETE SET NULL,
        sex              TEXT,
        dob              TEXT,
        genotype         TEXT,
        status           TEXT,
        notes            TEXT,
        mother_id        TEXT,
        father_id        TEXT,
        weight           NUMERIC,
        health_status    TEXT,
        strain           TEXT,
        source           TEXT,
        lab_id           TEXT,
        experiment_group TEXT,
        updated_at       TIMESTAMPTZ DEFAULT now()
      );")

    DBI::dbExecute(conn, "
      CREATE TABLE IF NOT EXISTS genes (
        gene     TEXT PRIMARY KEY,
        symbols  TEXT
      );")

    DBI::dbExecute(conn, "
      CREATE TABLE IF NOT EXISTS audit_log (
        id       BIGSERIAL PRIMARY KEY,
        time     TIMESTAMPTZ,
        action   TEXT,
        cage     TEXT,
        mouse_id TEXT,
        count    INT,
        details  TEXT
      );")
  }

  safe_read <- function(conn, tbl) {
    if (is.null(conn) || !tbl %in% DBI::dbListTables(conn)) return(NULL)
    DBI::dbReadTable(conn, tbl)
  }

  # Write-all strategy (simple + safe): replace table contents in a transaction.
  # Good for small/medium labs. For very large datasets, switch to UPSERTs.
  persist_all <- function(conn, cages, mice, genes, audit_append = NULL) {
    if (is.null(conn)) return()
    
    DBI::dbWithTransaction(conn, {
      if (!is.null(cages)) {
        DBI::dbExecute(conn, "DELETE FROM cages;")
        if (nrow(cages)) DBI::dbAppendTable(conn, "cages", cages)
      }
      if (!is.null(mice)) {
        DBI::dbExecute(conn, "DELETE FROM mice;")
        if (nrow(mice)) DBI::dbAppendTable(conn, "mice", mice)
      }
      if (!is.null(genes)) {
        DBI::dbExecute(conn, "DELETE FROM genes;")
        if (nrow(genes)) DBI::dbAppendTable(conn, "genes", genes)
      }
      if (!is.null(audit_append) && nrow(audit_append)) {
        DBI::dbAppendTable(conn, "audit_log", audit_append)
      }
    })
  }

  # One-time on startup: ensure schema, then load data into the reactiveVals
  init_db(pool)

  load_state_from_db <- function() {
    if (is.null(pool)) return(list(cages = NULL, mice = NULL, genes = NULL))
    
    list(
      cages = safe_read(pool, "cages"),
      mice  = safe_read(pool, "mice"),
      genes = safe_read(pool, "genes")
    )
  }

  # On stop: close pool
  onStop(function() {
    if (!is.null(pool)) pool::poolClose(pool)
  })

  # Enable a versioned S3 board (optional)
  s3_board <- tryCatch({
    if (nzchar(Sys.getenv("S3_BUCKET"))) {
      pins::board_s3(bucket = Sys.getenv("S3_BUCKET"),
                     region = Sys.getenv("S3_REGION"),
                     versioned = TRUE)
    } else {
      NULL
    }
  }, error = function(e) {
    message("S3 backup not configured: ", e$message)
    NULL
  })
  
  # --- helpers ---
  add_or_update_cage <- function(cage_name, genotype = "Unknown", sex = "Mixed",
                                 status = "Active", parent = "", setup = "") {
    cages <- cage_data_reactive()
    if (cage_name %in% cages$Cage_Name) {
      # update minimal fields if blank
      i <- which(cages$Cage_Name == cage_name)[1]
      if (nzchar(genotype) && cages$Genotype[i] == "") cages$Genotype[i] <- genotype
      if (nzchar(sex)      && cages$Sex[i]      == "") cages$Sex[i]      <- sex
      if (nzchar(status)   && cages$Status[i]   == "") cages$Status[i]   <- status
      if (nzchar(parent)   && cages$Parent_Cage[i] == "") cages$Parent_Cage[i] <- parent
      if (nzchar(setup)    && cages$Setup_Date[i] == "") cages$Setup_Date[i] <- setup
    } else {
      cages <- rbind(
        cages,
        data.frame(
          Cage_Name = cage_name,
          Genotype = genotype, Sex = sex, Num_Mice = 0, Num_Pups = 0,
          Status = status, Parent_Cage = parent,
          Genotyping_Done = FALSE, Genotype_Confirmed = FALSE,
          Setup_Date = setup, Notes = "", stringsAsFactors = FALSE
        )
      )
    }
    cage_data_reactive(cages)
  }
  
  recalc_cage_counts <- function() {
    cages <- cage_data_reactive()
    mice  <- mice_data_reactive()

    # total mice per cage
    total <- as.data.frame(table(mice$Cage), stringsAsFactors = FALSE)
    names(total) <- c("Cage_Name","Num_Mice_new")

    # pup mice per cage
    pup_df <- mice[is_pup_row(mice), , drop = FALSE]
    pup <- if (nrow(pup_df)) as.data.frame(table(pup_df$Cage), stringsAsFactors = FALSE)
           else data.frame(Var1 = character(0), Freq = integer(0), stringsAsFactors = FALSE)
    names(pup) <- c("Cage_Name","Num_Pups_new")

    cages <- cages %>%
      dplyr::left_join(total, by = "Cage_Name") %>%
      dplyr::left_join(pup,   by = "Cage_Name")

    cages$Num_Mice <- ifelse(is.na(cages$Num_Mice_new), 0L, cages$Num_Mice_new)
    cages$Num_Pups <- ifelse(is.na(cages$Num_Pups_new), 0L, cages$Num_Pups_new)
    cages$Num_Mice_new <- NULL
    cages$Num_Pups_new <- NULL

    cage_data_reactive(cages)
  }
  
  make_unique_ids <- function(prefix, n, existing) {
    out <- character(0); i <- 1
    while (length(out) < n) {
      cand <- sprintf("%s%02d", prefix, i)
      if (!(cand %in% existing)) out <- c(out, cand)
      i <- i + 1
    }
    out
  }
  
  distribute_round_robin <- function(vec, buckets) {
    if (!length(vec)) return(rep(NA_character_, 0))
    buckets[(seq_along(vec) - 1) %% length(buckets) + 1]
  }
  
  # Define %||% function (from rlang package)
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  # What counts as a pup (strict: only Status == "Pup")
  is_pup_row <- function(df) {
    s1 <- if ("Status" %in% names(df)) tolower(as.character(df$Status)) else ""
    s1 %in% "pup"
  }
  
  # Helper to count pups in a cage
  count_pups_in_cage <- function(cage_name) {
    m <- mice_data_reactive()
    m <- m[m$Cage == cage_name, , drop = FALSE]
    if (!nrow(m)) return(0)
    sum(is_pup_row(m), na.rm = TRUE)
  }
  
  # Helper to set cage sex from its mice
  set_cage_sex_from_mice <- function(cage_names) {
    cages <- cage_data_reactive()
    mice  <- mice_data_reactive()
    for (cn in cage_names) {
      sx <- unique(na.omit(as.character(mice$Sex[mice$Cage == cn])))
      new_sex <- if (length(sx) == 0) cages$Sex[cages$Cage_Name == cn]
                 else if (length(sx) == 1) sx
                 else "Mixed"
      cages$Sex[cages$Cage_Name == cn] <- ifelse(new_sex %in% c("M","F"), new_sex, "Mixed")
    }
    cage_data_reactive(cages)
  }
  
  # Gene catalog helpers
  split_syms <- function(s) gsub("^\\s+|\\s+$", "", unlist(strsplit(as.character(s %||% ""), ",")))
  
  # Normalize gene symbols (dedupe + trim)
  normalize_syms <- function(s) {
    paste(unique(trimws(unlist(strsplit(s %||% "", ",")))), collapse = ",")
  }
  get_gene_syms <- function(gene) {
    gc <- gene_catalog()
    if (!nrow(gc)) return(character(0))
    row <- gc[gc$Gene == gene, , drop = FALSE]
    if (!nrow(row)) character(0) else split_syms(row$Symbols[1])
  }
  
  # Update "Genotype" string by inserting/replacing one Gene(symbol) token
  # If s contains "Gene(...)" -> replace; else append.
  set_gene_in_string <- function(s, gene, sym) {
    s <- s %||% ""
    safe_gene <- gsub("([\\W])","\\\\\\1", gene, perl=TRUE)
    pat <- paste0("(?i)\\b", safe_gene, "\\s*\\([^)]*\\)")
    repl <- paste0(gene, "(", sym, ")")
    if (grepl(pat, s, perl = TRUE)) {
      gsub(pat, repl, s, perl = TRUE)
    } else if (nzchar(s)) {
      paste(s, repl, sep = " ; ")
    } else {
      repl
    }
  }
  
  # ---- Genotype parsing helpers (Gene(sym) in a single "Genotype" field) ----
  extract_gene_symbol <- function(genostr, gene) {
    if (is.null(genostr) || !nzchar(gene)) return(NA_character_)
    # find "Gene(...)"
    pat <- paste0("(?i)\\b", gsub("([\\W])","\\\\\\1", gene, perl=TRUE), "\\s*\\(([^)]*)\\)")
    m <- regexpr(pat, genostr, perl = TRUE)
    if (m[1] < 0) return(NA_character_)
    sub(pat, "\\1", regmatches(genostr, m))
  }
  
  has_gene_token <- function(genostr, gene) {
    !is.na(extract_gene_symbol(genostr, gene))
  }
  
  # Robust-ish DOB parser:
  # - YYYY-MM-DD -> Date
  # - M/D or M/D/YY(YY) -> assumes current year when year missing
  toDate_vec <- function(x) {
    x <- as.character(x %||% "")
    out <- rep(as.Date(NA), length(x))
    m1 <- grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
    m2 <- grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", x)
    m3 <- grepl("^\\d{1,2}/\\d{1,2}/\\d{2}$",  x)
    m4 <- grepl("^\\d{1,2}/\\d{1,2}$",         x)
    suppressWarnings({
      out[m1] <- as.Date(x[m1])
      out[m2] <- as.Date(x[m2], "%m/%d/%Y")
      out[m3] <- as.Date(x[m3], "%m/%d/%y")
      out[m4] <- as.Date(paste0(x[m4], "/", format(Sys.Date(), "%Y")), "%m/%d/%Y")
    })
    out
  }
  
  # Age calculator for Labmate export
  age_weeks <- function(dob_chr) {
    d <- toDate_vec(dob_chr)
    out <- floor(as.numeric(difftime(Sys.Date(), d, units = "days"))/7)
    ifelse(is.finite(out), out, NA_integer_)
  }
  
  # Shape a cage data.frame to your cage template columns (lab-mate style)
  cages_to_template <- function(df, cols) {
    out <- as.data.frame(
      setNames(replicate(length(cols), rep("", nrow(df)), simplify = FALSE), cols),
      stringsAsFactors = FALSE
    )
    common <- intersect(cols, names(df))
    for (nm in common) out[[nm]] <- as.character(df[[nm]])
    out
  }
  
  # Parse pasted table text (Excel -> TSV; CSV also works)
  read_pasted_df <- function(txt) {
    if (is.null(txt)) return(NULL)
    txt <- gsub("\r\n", "\n", txt, fixed = TRUE)
    txt <- gsub("\r", "\n", txt, fixed = TRUE)

    # Heuristic: Excel copy is usually tab-delimited
    sep <- if (grepl("\t", txt)) "\t" else ","

    con <- textConnection(txt)
    on.exit(close(con), add = TRUE)

    df <- tryCatch(
      utils::read.table(con, header = TRUE, sep = sep, quote = "\"",
                        check.names = FALSE, stringsAsFactors = FALSE,
                        comment.char = "", na.strings = c("", "NA")),
      error = function(e) NULL
    )
    if (is.null(df)) return(NULL)
    as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  # Smarter Excel reader (sheet/skip/header)
  read_upload_df <- function(path, sheet = NULL, skip = 0, col_names = TRUE) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("xlsx","xls")) {
      as.data.frame(
        readxl::read_excel(path, sheet = sheet %||% 1,
                           skip = skip %||% 0,
                           col_names = isTRUE(col_names)),
        stringsAsFactors = FALSE, check.names = FALSE
      )
    } else {
      utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
  }
  
  # Try to read only headers from mLIMS template once
  observe({
    path <- "mLIMS_Template.csv"           # <- place your template here
    if (file.exists(path)) {
      hdr <- names(utils::read.csv(path, nrows = 0, check.names = FALSE))
      if (length(hdr)) mlims_cols(hdr)
    } else {
      mlims_cols(NULL)  # fallback: no template found
    }
  })
  
  # Minimal, safe mapper: copy over only the columns you actually have.
  # Any template columns not present in your mice table are left blank.
  mice_to_mlims <- function(df, cols) {
    if (is.null(cols) || !length(cols)) return(df)           # fallback if no template
    # prepare empty data.frame with template columns
    out <- as.data.frame(
      setNames(
        replicate(length(cols), rep("", nrow(df)), simplify = FALSE),
        cols
      ),
      stringsAsFactors = FALSE
    )
    # copy overlapping columns verbatim (no fabrication)
    common <- intersect(cols, names(df))
    for (nm in common) out[[nm]] <- as.character(df[[nm]])
    out
  }
  
  # Dynamic Cage View UI renderer
  output$cage_view_ui <- renderUI({
    if (is.null(current_cage())) {
      # --- LIST STATE: show all cages ---
      tagList(
        fluidRow(
          column(12,
            div(class = "breadcrumb",
              span("🔎 Cage View", class = "breadcrumb-item")
            )
          )
        ),
        fluidRow(
          column(8,
            div(class = "search-container",
              textInput("cage_search_cv", "Search Cages (Press '/' to focus)",
                        placeholder = "Type to search cages...")
            )
          ),
          column(4,
            div(style = "margin-top:26px;",
              actionButton("reset_filters_cv", "Reset filters", class = "action-btn btn-edit")
            )
          )
        ),
        fluidRow(
          column(12,
            div(class="box", style="padding:12px; border:1px solid #eee; border-radius:8px; background:#fafafa;",
              h4("Filters"),
              fluidRow(
                column(3, uiOutput("cagecv_filter_gene_ui")),
                column(3, uiOutput("cagecv_filter_symbol_ui")),
                column(3, selectInput("cagecv_filter_sex", "Sex", c("Any","M","F","Mixed"), "Any")),
                column(3, selectInput("cagecv_filter_status", "Status",
                                      c("Any","Breeding","Holding","Experiment","Active","Retired"), "Any"))
              ),
              fluidRow(
                column(3, selectInput("cagecv_filter_has_pups", "Has pups", c("Any","Yes","No"), "Any")),
                column(3, selectInput("cagecv_filter_genotyped", "Genotyping done", c("Any","Yes","No"), "Any")),
                column(3, selectInput("cagecv_filter_confirmed", "Genotype confirmed", c("Any","Yes","No"), "Any")),
                column(3, dateRangeInput("cagecv_filter_setup", "Setup date", start = Sys.Date()-60, end = Sys.Date()))
              ),
              div(style="display:flex; gap:8px; margin-top:8px;",
                actionButton("cagecv_filters_reset", "Reset filters", class="action-btn btn-edit"),
                downloadButton("dl_filtered_cages_cv", "📥 Export filtered", class="action-btn btn-download"),
                div(style="flex:1"),
                checkboxInput("filter_genotyping_cv", "Show only cages with pending genotyping", value = FALSE)
              )
            )
          )
        ),
        fluidRow(
          column(12,
            div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;",
              span("💡 Tip: Click a row to open that cage. Use / to focus search.", class = "shortcut-hint"),
              div(
                actionButton("new_cage_from_cv", "➕ New Cage", class = "action-btn btn-view",
                           onclick = "Shiny.setInputValue('shortcut_new_cage', 1, {priority: 'event'})"),
                downloadButton("dl_all_cv", "📥 Download All Data", class = "action-btn btn-download"),
                actionButton("print_all_cv", "🖨️ Print All", class = "action-btn btn-download",
                           onclick = "Shiny.setInputValue('shortcut_print_all', 1, {priority: 'event'})"),
                actionButton("delete_selected_cages_cv", "🗑️ Delete Selected", class = "action-btn btn-danger")
              )
            )
          )
        ),
        fluidRow(
          column(12, DTOutput("cage_table_in_cage_view"))
        )
      )
    } else {
      # --- DETAIL STATE: show selected cage ---
      tagList(
        fluidRow(
          column(12,
            div(class = "breadcrumb",
              span("🔎 Cage View", class = "breadcrumb-item",
                   onclick = "Shiny.setInputValue('back_to_cage_list', 1, {priority:'event'})"),
              span(" > ", class = "breadcrumb-separator"),
              span(id = "current_cage_breadcrumb", paste0("Cage ", current_cage()), class = "breadcrumb-item")
            )
          )
        ),
        fluidRow(
          column(12,
            div(style = "display:flex; justify-content:space-between; align-items:center; margin-bottom:20px;",
              div(
                h2(textOutput("cage_view_title")),
                div(style = "font-size:14px; color:#666;",
                  span(textOutput("cage_header_genotype"), style="display:inline-block; margin-right:8px;"),
                  span(textOutput("cage_header_status"),   style="display:inline-block;")
                )
              ),
              div(
                actionButton("back_to_list_btn", "⬅️ Back to All Cages", class = "action-btn btn-edit"),
                actionButton("edit_cage", "✏️ Edit Cage", class = "action-btn btn-edit"),
                actionButton("add_mouse_to_cage", "➕ Add Mouse", class = "action-btn btn-view",
                             onclick = "Shiny.setInputValue('shortcut_add_mouse', 1, {priority: 'event'})"),
                downloadButton("dl_cage", "📥 Export Cage", class = "action-btn btn-download"),
                downloadButton("dl_cage_csv", "⬇️ CSV (mLIMS)", class = "action-btn btn-download"),
                downloadButton("dl_labmate_cage", "📥 Labmate (this cage)", class = "action-btn btn-download"),
                actionButton("print_current_cage", "🖨️ Print Cage", class = "action-btn btn-download",
                             onclick = "Shiny.setInputValue('shortcut_print_cage', 1, {priority: 'event'})"),
                actionButton("delete_cage", "🗑️ Delete This Cage", class = "action-btn btn-danger")
              )
            )
          )
        ),
        fluidRow(
          column(12,
            div(style = "margin-bottom:20px; padding:15px; background-color:#f8f9fa; border-radius:8px;",
              h4("Quick Actions"),
              div(style = "display:flex; gap:10px; flex-wrap:wrap;",
                actionButton("add_pup_to_cage", "🐭 Add Pup", class = "action-btn btn-view"),
                actionButton("mark_weaned", "📦 Wean / Split", class = "action-btn btn-edit"),
                actionButton("split_cage", "🔀 Split Cage", class = "action-btn btn-edit"),
                actionButton("add_genotype", "🧬 Add Genotype", class = "action-btn btn-download"),
                actionButton("move_cage", "🏠 Move to New Cage", class = "action-btn btn-edit")
              )
            )
          )
        ),
        
        # Template control bar
        fluidRow(
          column(12,
            div(style = "margin: 10px 0 0 0; padding: 10px; background:#f8f9fa; border-radius:8px;",
              checkboxInput("use_mlims_template", "Show as mLIMS CSV template", value = TRUE),
              span("If the template file isn't found, the app uses the mice table as-is.", style="color:#777; font-size:12px; margin-left:6px;")
            )
          )
        ),
        fluidRow(
          column(12,
            div(style = "margin-bottom:20px; padding:15px; background-color:#f8f9fa; border-radius:8px;",
              h4("Cage Information"),
              verbatimTextOutput("cage_metadata")
            )
          )
        ),
        fluidRow(
          column(12,
            div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
              h4("Mice in this Cage"),
              actionButton("delete_selected_mice_in_cage", "🗑️ Delete Selected", class = "action-btn btn-danger", style="margin-left:10px;")
            ),
            conditionalPanel("input.use_mlims_template", DTOutput("mice_in_cage_table_mlims")),
            conditionalPanel("!input.use_mlims_template", DTOutput("mice_in_cage_table"))
          )
        ),
        fluidRow(
          column(12,
            h4("Raw CSV Preview"),
            tags$pre(style="background:#fff;border:1px solid #eee;padding:10px;overflow:auto;",
                     textOutput("cage_csv_preview"))
          )
        )
      )
    }
  })
  
  # Upload modal UI helper function
  uploadModalUI <- function(has_template) {
    modalDialog(
      title = "Upload / Import mLIMS Data",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close"),
        actionButton("commit_import", "✅ Commit Import", class = "action-btn btn-view", id = "commit_import_btn")
      ),
      div(
        # What are we importing?
        radioButtons("import_kind", "What are you importing?",
                     choices = c("Mouse-based (one row per animal)" = "mice",
                                 "Cage-based (one row per cage)" = "cage"),
                     inline = TRUE),

        # SOURCE: file upload OR paste table
        tabsetPanel(id = "import_source",
          tabPanel("Upload file",
            br(),
            fileInput("upload_file", "Select CSV or Excel (XLSX/XLS) file",
                      accept = c(".csv", ".xlsx", ".xls")),
            uiOutput("excel_options_ui"),  # shows only for Excel files
            div(style="display:flex; gap:16px; align-items:center; margin-top:8px;",
                sliderInput("preview_rows", "Max preview rows", min = 10, max = 2000, value = 200, step = 10, width = "380px"),
                tags$small(style="color:#777;", "Large files: preview is truncated to speed things up.")
            ),
            tags$small("Tip: Excel is convenient; CSV is fastest for multi-hundred-thousand row imports.")
          ),
          tabPanel("Paste table",
            br(),
            tags$small("From Excel or Google Sheets: select your table (include the header row), press Ctrl/Cmd+C, then paste here."),
            textAreaInput("paste_table", NULL, placeholder = "Paste your table here (tab- or comma-delimited)", rows = 10),
            actionLink("clear_paste", "Clear")
          )
        ),

        # Template downloads
        div(style="display:flex; gap:8px; align-items:center; margin:12px 0;",
          downloadButton("dl_mlims_mouse_template", "⬇️ Mouse Template (CSV)", class = "action-btn btn-download"),
          downloadButton("dl_mlims_cage_template",  "⬇️ Cage Template (CSV)",  class = "action-btn btn-download"),
          span(if (has_template) "Mouse template uses your mLIMS headers." else
                 "Place mLIMS_Template.csv next to the app to use your headers.",
               style="color:#777; font-size:12px; margin-left:6px;")
        ),

        # Options
        checkboxInput("import_upsert", "Update existing records with uploaded (non-empty) values", TRUE),
        checkboxInput("import_recount", "Recount Num_Mice in cages from imported mice", TRUE),

        hr(),
        h4("Preview"),
        DT::DTOutput("import_preview"),
        hr(),
        h4("Issues"),
        htmlOutput("import_issues")
      )
    )
  }
  
  # Cage overview table with checkboxes and bulk delete
  output$cage_overview_table <- DT::renderDT({
    df <- filtered_cages_dashboard()

    # Add/keep Has_Pups and column order for display
    df$Has_Pups <- ifelse(df$Num_Pups > 0, "Yes", "No")
    df <- df[, c("Cage_Name","Genotype","Sex","Num_Mice","Num_Pups","Has_Pups",
                 "Status","Parent_Cage","Genotyping_Done","Genotype_Confirmed","Setup_Date"),
             drop = FALSE]

    # clickable parent link
    df$Parent_Cage <- ifelse(
      nzchar(df$Parent_Cage),
      sprintf("<a href='#' class='parent-cage' data-cage='%s'>%s</a>", df$Parent_Cage, df$Parent_Cage),
      ""
    )

    # add leading blank col for select checkbox
    df2 <- cbind(` ` = "", df)

    datatable(
      df2,
      extensions = c("Buttons","Select"),
      selection = "multiple",
      escape = FALSE, rownames = FALSE,
      options = list(
        dom = 'Bfrtip',
        buttons = c('copy','csv','excel','pdf','print'),
        pageLength = 10,
        order = list(list(1,'asc')), # sort by Cage_Name (now at col 1)
        columnDefs = list(
          list(orderable = FALSE, className = 'select-checkbox', targets = 0)
        ),
        select = list(style = 'multi', selector = 'td:first-child'),
        server = FALSE
      ),
      callback = JS("
        // open parent cage
        table.on('click', 'a.parent-cage', function(e){
          e.preventDefault();
          var cage = $(this).data('cage');
          Shiny.setInputValue('parent_cage_click', cage, {priority:'event'});
        });
        // open cage when clicking its name (col 1)
        table.on('click', 'td:nth-child(2)', function(){
          var data = table.row(this).data();
          Shiny.setInputValue('cage_click', data[1], {priority:'event'});
        });
        // send selected Cage_Name keys to Shiny
        table.on('select.dt deselect.dt', function(e, dt){
          var cages = dt.rows({selected:true}).data().toArray().map(function(r){ return r[1]; });
          Shiny.setInputValue('cage_selected_keys_dashboard', cages, {priority:'event'});
        });
      ")
    ) %>%
      formatStyle('Cage_Name', cursor='pointer', color='#2196f3', fontWeight='bold') %>%
      formatStyle('Parent_Cage', cursor='pointer', color='#ff5722', fontWeight='bold') %>%
      formatStyle('Status',
        backgroundColor = styleEqual(
          c('Breeding','Holding','Experiment','Active','Retired'),
          c('#ffeb3b','#2196f3','#ff9800','#4caf50','#9e9e9e')
        )
      ) %>%
      formatStyle('Genotyping_Done', backgroundColor = styleEqual(c(TRUE,FALSE), c('#4caf50','#f44336'))) %>%
      formatStyle('Genotype_Confirmed', backgroundColor = styleEqual(c(TRUE,FALSE), c('#4caf50','#f44336'))) %>%
      formatStyle('Has_Pups',
        backgroundColor = styleEqual(c('Yes','No'), c('#e8f5e9','#ffebee')),
        color = styleEqual(c('Yes','No'), c('#1b5e20','#b71c1c')),
        fontWeight = 'bold'
      )
  })
  
  # Reactive header outputs for Cage View detail state
  output$cage_view_title <- renderText({
    if (is.null(current_cage())) "All Cages" else paste("Cage:", current_cage())
  })
  
  output$cage_header_genotype <- renderText({
    req(current_cage())
    info <- cage_data_reactive()[cage_data_reactive()$Cage_Name == current_cage(), ]
    if (nrow(info)) paste0("Genotype: ", info$Genotype[1]) else "Genotype: "
  })
  
  output$cage_header_status <- renderText({
    req(current_cage())
    info <- cage_data_reactive()[cage_data_reactive()$Cage_Name == current_cage(), ]
    if (nrow(info)) paste0("| Status: ", info$Status[1]) else "| Status: "
  })
  
  # Cage metadata display
  output$cage_metadata <- renderPrint({
    if (is.null(current_cage())) {
      cat("No cage selected. Use the table to pick one.")
      return(invisible(NULL))
    }
    cage_info <- cage_data_reactive()[cage_data_reactive()$Cage_Name == current_cage(), ]
    if (nrow(cage_info) > 0) {
      cat("Cage Name:", cage_info$Cage_Name[1], "\n")
      cat("Genotype:", cage_info$Genotype[1], "\n")
      cat("Sex:", cage_info$Sex[1], "\n")
      cat("Number of Mice:", cage_info$Num_Mice[1], "\n")
      cat("Pup Count:", cage_info$Num_Pups[1], "\n")
      cat("Has Pups:", ifelse(cage_info$Num_Pups[1] > 0, "Yes", "No"), "\n")
      cat("Status:", cage_info$Status[1], "\n")
      if (nzchar(cage_info$Parent_Cage[1])) cat("Parent Cage:", cage_info$Parent_Cage[1], "\n")
      cat("Genotyping Done:", cage_info$Genotyping_Done[1], "\n")
      cat("Genotype Confirmed:", cage_info$Genotype_Confirmed[1], "\n")
      if (nzchar(cage_info$Setup_Date[1])) cat("Setup Date:", cage_info$Setup_Date[1], "\n")
    }
  })
  
  # Mice in cage table with checkboxes and bulk delete
  output$mice_in_cage_table_mlims <- renderDT({
    # pick mice for the current cage (or show a friendly message)
    if (is.null(current_cage())) {
      return(datatable(
        data.frame(Message = "Select a cage to see its mice."),
        options = list(dom = 't'),
        rownames = FALSE
      ))
    }

    df <- mice_data_reactive()[mice_data_reactive()$Cage == current_cage(), , drop = FALSE]

    # If the user wants template view AND we have a template, reshape to its columns.
    if (isTRUE(input$use_mlims_template) && !is.null(mlims_cols())) {
      df <- mice_to_mlims(df, mlims_cols())
    }

    # add leading blank col for select checkbox
    df2 <- cbind(` ` = "", df)

    datatable(
      df2,
      extensions = c("Buttons","Select"),
      selection = "multiple", escape = FALSE, rownames = FALSE,
      options = list(
        pageLength = 10,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
        columnDefs = list(list(orderable=FALSE, className='select-checkbox', targets=0)),
        select = list(style='multi', selector='td:first-child')
      ),
      callback = JS("
        // whenever selection changes, send selected Mouse_IDs
        table.on('select.dt deselect.dt', function(e, dt){
          // col 2 is Mouse_ID because col 1 is the checkbox
          var ids = dt.rows({selected:true}).data().toArray().map(function(r){ return r[1]; });
          Shiny.setInputValue('selected_mouse_ids', ids, {priority:'event'});
        });
      ")
    )
  }, server = FALSE)
  
  # Regular mice in cage table (for weaning selection)
  output$mice_in_cage_table <- DT::renderDT({
    req(current_cage())
    cage_mice <- mice_data_reactive()[mice_data_reactive()$Cage == current_cage(), ]
    df2 <- cbind(` ` = "", cage_mice)
    datatable(
      df2,
      extensions = c("Buttons","Select"),
      selection = "multiple", escape = FALSE, rownames = FALSE,
      options = list(
        dom = 'Bfrtip',
        buttons = c('copy','csv','excel','pdf','print'),
        pageLength = 10,
        columnDefs = list(list(orderable = FALSE, className = 'select-checkbox', targets = 0)),
        select = list(style = 'multi', selector = 'td:first-child')
      ),
      callback = JS("
        // open detail when clicking Mouse_ID (2nd visible column)
        table.on('click', 'td:nth-child(2)', function(){
          var data = table.row(this).data();
          if (data && data.length > 1) {
            Shiny.setInputValue('mouse_click', data[1], {priority:'event'});
          }
        });

        // whenever selection changes, send selected Mouse_IDs
        table.on('select.dt deselect.dt', function(e, dt){
          // col 2 is Mouse_ID because col 1 is the checkbox
          var ids = dt.rows({selected:true}).data().toArray().map(function(r){ return r[1]; });
          Shiny.setInputValue('selected_mouse_ids', ids, {priority:'event'});
        });
      ")
    ) %>%
      formatStyle('Mouse_ID', cursor='pointer', color='#2196f3', fontWeight='bold')
  }, server = FALSE)
  
  # Mouse table with row selection
  # MAIN mice table with checkboxes + bulk delete
  output$mice_table <- DT::renderDT({
    df <- filtered_mice()
    df2 <- cbind(` ` = "", df)  # checkbox column

    datatable(
      df2,
      extensions = c("Buttons","Select"),
      selection = "multiple", escape = FALSE, rownames = FALSE,
      options = list(
        dom = 'Bfrtip', buttons = c('copy','csv','excel','pdf','print'),
        pageLength = 10, order = list(list(2,'asc')), # Mouse_ID col now #2
        columnDefs = list(list(orderable=FALSE, className='select-checkbox', targets=0)),
        select = list(style='multi', selector='td:first-child'),
        server = FALSE
      ),
      callback = JS("
        // open detail by clicking Mouse_ID (col 2)
        table.on('click', 'td:nth-child(2)', function(){
          var data = table.row(this).data();
          if (data && data.length > 1) {
            Shiny.setInputValue('mouse_click', data[1], {priority:'event'});
          }
        });
        // send selected Mouse_ID keys to Shiny
        table.on('select.dt deselect.dt', function(e, dt){
          var ids = dt.rows({selected:true}).data().toArray().map(function(r){ return r[1]; });
          Shiny.setInputValue('mouse_selected_keys_main', ids, {priority:'event'});
        });
      ")
    ) %>%
      formatStyle('Mouse_ID', cursor='pointer', color='#2196f3', fontWeight='bold')
  })
  
  # Handle cage table row selection - REMOVED to allow checkbox selection for bulk operations
  # Navigation is handled by cage_click callback instead
  
  # Handle cage clicks - navigate to cage view (for backward compatibility)
  observeEvent(input$cage_click, {
    req(input$cage_click)
    current_cage(input$cage_click)
    current_mouse(NULL)
    updateTabItems(session, "sidebar", "cage_view")
  })
  
  # Handle mouse table row selection (from main mouse table)
  observeEvent(input$mice_table_rows_selected, {
    # We rely on clicking Mouse_ID (callback) instead of row selection now.
  })
  
  # Handle mouse clicks from main table
  observeEvent(input$mouse_click, {
    current_mouse(input$mouse_click)
    mr <- mice_data_reactive()[mice_data_reactive()$Mouse_ID == input$mouse_click, , drop = FALSE]
    showMouseDetailModal(mr)
  })
  

  
  # Function to show mouse detail modal
  showMouseDetailModal <- function(mouse_row) {
    if (is.null(mouse_row) || !nrow(mouse_row)) return(invisible(NULL))
    mr <- mouse_row[1, , drop = FALSE]
    mid <- as.character(mr$Mouse_ID)
    mom <- as.character(mr$Mother_ID %||% "")
    dad <- as.character(mr$Father_ID %||% "")
    cag <- as.character(mr$Cage)

    link_or_text <- function(id, label) {
      if (!nzchar(id)) return(span("—"))
      tags$a(href = "#", onclick = sprintf("Shiny.setInputValue('mouse_click','%s',{priority:'event'})", id), label)
    }

    showModal(modalDialog(
      title = paste("Mouse:", mid),
      easyClose = TRUE,
      size = "m",
      div(
        div(class = "breadcrumb",
            span("🐭 Mouse Records", class = "breadcrumb-item",
                 onclick = "Shiny.setInputValue('breadcrumb_click', 'mice', {priority:'event'})"),
            span(" > ", class = "breadcrumb-separator"),
            span(paste("Mouse", mid), class = "breadcrumb-item")
        ),
        tags$div(style="display:grid; grid-template-columns: 160px 1fr; gap:6px; align-items:center;",
          tags$strong("Mouse ID:"), span(mid),
          tags$strong("Cage:"), tags$a(href="#",
            onclick = sprintf("Shiny.setInputValue('cage_click','%s',{priority:'event'})", cag), cag),
          tags$strong("Sex:"), span(mr$Sex %||% ""),
          tags$strong("DOB:"), span(mr$DOB %||% ""),
          tags$strong("Genotype:"), span(mr$Genotype %||% ""),
          tags$strong("Status:"), span(mr$Status %||% ""),
          tags$strong("Mother:"), link_or_text(mom, mom),
          tags$strong("Father:"), link_or_text(dad, dad),
          tags$strong("Strain:"), span(mr$Strain %||% ""),
          tags$strong("Source:"), span(mr$Source %||% ""),
          tags$strong("Weight (g):"), span(mr$Weight %||% ""),
          tags$strong("Health:"), span(mr$Health_Status %||% ""),
          tags$strong("Lab ID:"), span(mr$Lab_ID %||% ""),
          tags$strong("Experiment Group:"), span(mr$Experiment_Group %||% ""),
          tags$strong("Notes:"), span(mr$Notes %||% "")
        ),
        br(),
        div(style="text-align:center;",
          actionButton("edit_mouse_modal", "✏️ Edit Mouse", class = "action-btn btn-edit"),
          actionButton("delete_mouse", "🗑️ Delete Mouse", class = "action-btn btn-danger")
        )
      )
    ))
  }
  
  # Helper to show a quick "Add Gene" modal
  showQuickAddGeneModal <- function() {
    showModal(modalDialog(
      title = "➕ Add a Gene",
      easyClose = TRUE,
      size = "m",
      div(
        textInput("quick_gene_name", "Gene name:", placeholder = "e.g., Cx3cr1"),
        textInput("quick_gene_syms", "Allowed symbols (comma-separated):",
                  placeholder = "e.g., +/+, +/-, -/-, WT, flox/+, flox/flox"),
        tags$small("Tip: You can use '/' freely now (e.g., flox/+).")
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("quick_add_gene_save", "Save & Continue", class = "action-btn btn-view")
      )
    ))
  }
  
  # Helper to show genotype modal with optional pre-selected gene
  showGenotypeModal <- function(pref_gene = NULL) {
    gc <- gene_catalog()
    has_genes <- nrow(gc) > 0
    default_gene <- pref_gene %||% last_gene() %||% if (has_genes) gc$Gene[1] else ""

    showModal(modalDialog(
      title = paste("Add Genotype —", if (is.null(current_mouse())) current_cage() else current_mouse()),
      size = "m",
      easyClose = TRUE,
      div(
        if (!has_genes) div(style="margin-bottom:8px; color:#c62828;",
                            "No genes yet. Please add one."),
        selectInput("apply_scope", "Apply to:",
                    choices = c("Selected mice in this cage" = "sel",
                                "All mice in this cage"      = "all",
                                "Cage record only"           = "cage"),
                    selected = "all"),
        selectizeInput("gene_choice", "Gene:",
                       choices = if (has_genes) gc$Gene else character(0),
                       selected = if (has_genes) default_gene else NULL,
                       options = list(create = FALSE)),
        uiOutput("gene_symbol_ui"),
        dateInput("genotype_date", "Genotyping Date:", value = Sys.Date()),
        textInput("genotype_lab", "Lab / Assay:", value = ""),
        textInput("genotype_notes", "Notes:", value = "")
      ),
      footer = tagList(
        actionButton("goto_genes_tab", "➕ Add Gene", class="action-btn btn-edit"),
        modalButton("Cancel"),
        actionButton("save_genotype_new", "Save", class = "action-btn btn-download",
                     disabled = !has_genes)
      )
    ))
  }
  
  # Reusable Add-Mouse modal: preselect a cage (if given), or let user pick/create
  showAddMouseModal <- function(default_cage = NULL) {
    cages <- sort(unique(cage_data_reactive()$Cage_Name))

    showModal(modalDialog(
      title = if (is.null(default_cage)) "Add Mouse" else paste("Add Mouse to", default_cage),
      easyClose = TRUE, size = "m",
      div(
        textInput("new_mouse_id_unified", "Mouse ID:"),
        selectizeInput(
          "new_mouse_cage_select",
          "Cage (select existing or type a new name):",
          choices  = cages,
          selected = default_cage,
          options  = list(create = TRUE, placeholder = "Choose a cage or type a new one")
        ),
        selectInput("new_mouse_sex_unified", "Sex:", choices = c("M","F")),
        textInput("new_mouse_dob_unified", "DOB (YYYY-MM-DD):"),
        textInput("new_mouse_genotype_unified", "Genotype:", value = "Unknown"),
        selectInput("new_mouse_status_unified", "Status:",
                    c("Pup","Juvenile","Active","Holding","Experiment","Retired"), selected = "Active"),
        numericInput("new_mouse_weight_unified", "Weight (g):", value = NA, min = 0),
        textInput("new_mouse_health_unified", "Health Status:", value = ""),
        textInput("new_mouse_strain_unified", "Strain:", value = ""),
        textInput("new_mouse_source_unified", "Source:", value = ""),
        textInput("new_mouse_labid_unified", "Lab ID:", value = ""),
        textInput("new_mouse_group_unified", "Experiment Group:", value = ""),
        textInput("new_mouse_notes_unified", "Notes:", value = "")
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_new_mouse_unified", "Save", class = "action-btn btn-view")
      )
    ))
  }
  

  
  # Handle breadcrumb navigation
  observeEvent(input$breadcrumb_click, {
    if (input$breadcrumb_click == "home") {
      updateTabItems(session, "sidebar", "dashboard")
      current_cage(NULL)
    }
  })
  
  # Keyboard shortcut handlers
  observeEvent(input$shortcut_new_cage, {
    showModal(modalDialog(
      title = "Add New Cage",
      div(
        textInput("new_cage_id", "Cage ID:"),
        selectInput("new_cage_status", "Status:",
                    choices = c("Breeding","Holding","Experiment","Active","Retired"),
                    selected = "Active"),
        selectInput("new_cage_genotype", "Genotype:", 
                   choices = c("Cx3cr1", "Tmem119", "Other")),
        selectInput("new_cage_sex", "Sex:", choices = c("M", "F", "Mixed")),
        textInput("new_cage_notes", "Notes:")
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_new_cage", "Save", class = "action-btn btn-view")
      )
    ))
  })
  
  # Download handlers
  output$dl_all <- downloadHandler(
    filename = function() sprintf("mLIMS_export_%s.xlsx", Sys.Date()),
    content = function(file) {
      # Mouse sheet: use mLIMS mouse headers if available
      mice_df <- mice_data_reactive()
      if (!is.null(mlims_cols())) mice_df <- mice_to_mlims(mice_df, mlims_cols())

      # Cage sheet: use your cage template columns
      cages_df <- cages_to_template(cage_data_reactive(), cage_template_cols)

      writexl::write_xlsx(list(
        Cages = cages_df,
        Mice  = mice_df
      ), path = file)
    }
  )
  
  output$dl_cage <- downloadHandler(
    filename = function() sprintf("mLIMS_cage_%s_%s.xlsx", current_cage(), Sys.Date()),
    content = function(file) {
      cage_info <- cage_data_reactive()[cage_data_reactive()$Cage_Name == current_cage(), ]
      cage_mice <- mice_data_reactive()[mice_data_reactive()$Cage == current_cage(), ]

      # Shape to templates
      cage_info_tmpl <- cages_to_template(cage_info, cage_template_cols)
      cage_mice_tmpl <- if (!is.null(mlims_cols())) mice_to_mlims(cage_mice, mlims_cols()) else cage_mice

      writexl::write_xlsx(list(
        Cage_Info = cage_info_tmpl,
        Mice      = cage_mice_tmpl
      ), path = file)
    }
  )
  
  # CSV preview output
  output$cage_csv_preview <- renderText({
    if (is.null(current_cage())) return("Select a cage to see CSV.")
    df <- mice_data_reactive()[mice_data_reactive()$Cage == current_cage(), , drop = FALSE]
    if (isTRUE(input$use_mlims_template) && !is.null(mlims_cols())) {
      df <- mice_to_mlims(df, mlims_cols())
    }
    tc <- textConnection("x", "w", local = TRUE)
    utils::write.csv(df, tc, row.names = FALSE, na = "")
    close(tc)
    paste(x, collapse = "\n")
  })
  
  # mLIMS CSV download handler
  output$dl_cage_csv <- downloadHandler(
    filename = function() sprintf("mLIMS_cage_%s_%s.csv", current_cage(), Sys.Date()),
    content = function(file) {
      df <- mice_data_reactive()[mice_data_reactive()$Cage == current_cage(), , drop = FALSE]
      if (isTRUE(input$use_mlims_template) && !is.null(mlims_cols())) {
        df <- mice_to_mlims(df, mlims_cols())
      }
      utils::write.csv(df, file, row.names = FALSE, na = "")
    }
  )
  
  # Download handler for Cage View (same as main)
  output$dl_all_cv <- downloadHandler(
    filename = function() sprintf("mLIMS_export_%s.xlsx", Sys.Date()),
    content = function(file) {
      # Mouse sheet: use mLIMS mouse headers if available
      mice_df <- mice_data_reactive()
      if (!is.null(mlims_cols())) mice_df <- mice_to_mlims(mice_df, mlims_cols())

      # Cage sheet: use your cage template columns
      cages_df <- cages_to_template(cage_data_reactive(), cage_template_cols)

      writexl::write_xlsx(list(
        Cages = cages_df,
        Mice  = mice_df
      ), path = file)
    }
  )
  
  # Template download handlers
  output$dl_mlims_mouse_template <- downloadHandler(
    filename = function() "mLIMS_Mouse_Template.csv",
    content = function(file) {
      cols <- if (!is.null(mlims_cols())) mlims_cols() else names(mice_data_reactive())
      empty <- as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols))
      utils::write.csv(empty, file, row.names = FALSE, na = "")
    }
  )
  
  output$dl_mlims_cage_template <- downloadHandler(
    filename = function() "mLIMS_Cage_Template.csv",
    content = function(file) {
      empty <- as.data.frame(setNames(replicate(length(cage_template_cols), character(0), simplify = FALSE),
                                      cage_template_cols))
      utils::write.csv(empty, file, row.names = FALSE, na = "")
    }
  )
  
  # All mice -> mLIMS-like Excel
  output$dl_mice_all_mlims_xlsx <- downloadHandler(
    filename = function() sprintf("mLIMS_mice_all_%s.xlsx", Sys.Date()),
    content = function(file) {
      df <- mice_data_reactive()
      if (!is.null(mlims_cols())) df <- mice_to_mlims(df, mlims_cols())
      writexl::write_xlsx(list(Mice = df), path = file)
    }
  )
  
  # Filtered mice -> mLIMS-like Excel
  output$dl_mice_filtered_mlims_xlsx <- downloadHandler(
    filename = function() sprintf("mLIMS_mice_filtered_%s.xlsx", Sys.Date()),
    content = function(file) {
      df <- filtered_mice()
      if (!is.null(mlims_cols())) df <- mice_to_mlims(df, mlims_cols())
      writexl::write_xlsx(list(Mice = df), path = file)
    }
  )
  
  # Cages-only Excel export
  output$dl_cages_mlims_xlsx <- downloadHandler(
    filename = function() sprintf("mLIMS_cages_%s.xlsx", Sys.Date()),
    content = function(file) {
      cages <- cage_data_reactive()
      mice  <- mice_data_reactive()
      # keep only cages that still have mice rows (avoids seeing a name you deleted earlier in another session)
      keep <- cages$Cage_Name %in% unique(mice$Cage)
      cages <- cages[keep, , drop = FALSE]
      df <- cages_to_template(cages, cage_template_cols)
      writexl::write_xlsx(list(Cages = df), path = file)
    }
  )
  
  # Labmate format exports
  output$dl_labmate_all <- downloadHandler(
    filename = function() sprintf("mice_labmate_%s.xlsx", Sys.Date()),
    content = function(file) {
      m <- mice_data_reactive()
      need <- c("Cage","Mouse_ID","DOB","Sex","Genotype")
      for (nm in setdiff(need, names(m))) m[[nm]] <- NA

      out <- data.frame(
        `Pack cages` = m$Cage,
        `Mouse#`     = m$Mouse_ID,
        `DOB`        = m$DOB,
        `Age`        = age_weeks(m$DOB),
        `Sex`        = m$Sex,
        `Genotype`   = m$Genotype,
        `DOB.1`      = "",          # left blank (not tracked)
        `Toe cutting`= "",          # left blank (not tracked)
        `Wean`       = "",          # left blank (not tracked)
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      writexl::write_xlsx(list(`Mice` = out), path = file)
    }
  )

  output$dl_labmate_cage <- downloadHandler(
    filename = function() sprintf("mice_labmate_%s_%s.xlsx", gsub("[^A-Za-z0-9]+","_", current_cage()), Sys.Date()),
    content = function(file) {
      req(current_cage())
      m <- mice_data_reactive()
      m <- m[m$Cage == current_cage(), , drop=FALSE]
      need <- c("Cage","Mouse_ID","DOB","Sex","Genotype")
      for (nm in setdiff(need, names(m))) m[[nm]] <- NA

      out <- data.frame(
        `Pack cages` = m$Cage,
        `Mouse#`     = m$Mouse_ID,
        `DOB`        = m$DOB,
        `Age`        = age_weeks(m$DOB),
        `Sex`        = m$Sex,
        `Genotype`   = m$Genotype,
        `DOB.1`      = "",
        `Toe cutting`= "",
        `Wean`       = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      writexl::write_xlsx(list(`Mice` = out), path = file)
    }
  )

  # Helper: only keep mice that belong to existing (non-deleted) cages
  live_mice <- reactive({
    m <- mice_data_reactive()
    live_cages <- unique(cage_data_reactive()$Cage_Name)
    m[m$Cage %in% live_cages, , drop = FALSE]
  })

  output$dl_labmate_xlsx <- downloadHandler(
    filename = function() sprintf("Labmate_mice_%s.xlsx", Sys.Date()),
    content = function(file) {
      m <- live_mice()

      # Robust DOB -> Date
      dob <- toDate_vec(m$DOB)
      age_days <- as.integer(difftime(Sys.Date(), dob, units = "days"))
      # Show days when close to birth; otherwise months (matches your 4 / 3 in the screenshot)
      age_display <- ifelse(
        is.na(age_days), NA,
        ifelse(abs(age_days) < 60, age_days, floor(age_days / 30))
      )

      lab <- data.frame(
        `Pack cages`  = m$Cage,
        `Mouse#`      = m$Mouse_ID,
        `DOB`         = m$DOB,
        `Age`         = age_display,
        `Sex`         = m$Sex,
        `Genotype`    = m$Genotype,
        `DOB.1`       = "",            # left blank to match your labmate sheet
        `Toe cutting` = "",
        `Wean`        = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      # Sort like the example: by cage then mouse id
      lab <- lab[order(lab$`Pack cages`, lab$`Mouse#`), ]

      writexl::write_xlsx(list(`Mice` = lab), path = file)
    }
  )
  
  # File reading and validation functions
  read_upload_df <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("xlsx","xls")) {
      as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
  }
  
  validate_mice_df <- function(df) {
    errs <- c(); warns <- c()
    # Key columns we can accept from your template; we won't fabricate if absent.
    key_candidates <- c("Mouse_ID","Animal_ID","Physical Tag","Earmark")
    key_present <- key_candidates[key_candidates %in% names(df)]
    if (!length(key_present)) errs <- c(errs, "Missing a mouse ID column (e.g., Mouse_ID / Animal_ID / Physical Tag / Earmark).")
    if (!"Cage" %in% names(df)) warns <- c(warns, "Column 'Cage' is missing; mice will not be linked to cages.")
    list(errors = errs, warnings = warns)
  }
  
  validate_cage_df <- function(df) {
    errs <- c(); warns <- c()
    if (!"Cage_Name" %in% names(df)) errs <- c(errs, "Missing required 'Cage_Name'.")
    list(errors = errs, warnings = warns)
  }
  
  # Cage selector functionality
  # When you navigate to Cage View, clear current cage to show the list
  observeEvent(input$sidebar, {
    if (identical(input$sidebar, "cage_view")) {
      current_cage(NULL)  # show the list first
    }
  })
  

  
  # Disable "Export Cage" when no cage is selected (nice UX)
  observe({
    shinyjs::toggleState("dl_cage", !is.null(current_cage()))
  })
  
  # Disable CSV button when no cage is selected
  observe({
    shinyjs::toggleState("dl_cage_csv", !is.null(current_cage()))
  })
  
  # Row click handler for Cage View list table - REMOVED to allow checkbox selection for bulk operations
  # Navigation is handled by cage_click callback instead
  
  # Reset filters button for Cage View
  observeEvent(input$reset_filters_cv, {
    updateTextInput(session, "cage_search_cv", value = "")
    updateCheckboxInput(session, "filter_genotyping_cv", value = FALSE)
  })
  
  # File upload processing
  observeEvent(input$upload_file, {
    req(input$upload_file$datapath)
    df <- read_upload_df(input$upload_file$datapath)

    if (identical(input$import_kind, "mice")) {
      # Shape to your mLIMS template header if available; leave unmapped columns blank.
      if (!is.null(mlims_cols())) {
        cols <- mlims_cols()
        out <- as.data.frame(
          setNames(replicate(length(cols), rep("", nrow(df)), simplify = FALSE), cols),
          stringsAsFactors = FALSE
        )
        common <- intersect(names(df), cols)
        for (nm in common) out[[nm]] <- as.character(df[[nm]])
        df <- out
      }
      v <- validate_mice_df(df)
      upload_issues(v)
      uploaded_df(df)
    } else {
      # cage-based upload
      # Keep only columns we know in the cage template (extra columns are ignored)
      keep <- intersect(names(df), cage_template_cols)
      df2 <- as.data.frame(df[keep], stringsAsFactors = FALSE, check.names = FALSE)
      # Add any missing template cols as blanks
      missing <- setdiff(cage_template_cols, names(df2))
      for (m in missing) df2[[m]] <- ""
      df2 <- df2[, cage_template_cols, drop = FALSE]
      v <- validate_cage_df(df2)
      upload_issues(v)
      uploaded_df(df2)
    }
  })
  
  # Cage View list table with checkboxes and bulk delete
  output$cage_table_in_cage_view <- DT::renderDT({
    df <- filtered_cages_cv()

    df$Has_Pups <- ifelse(df$Num_Pups > 0, "Yes", "No")
    df <- df[, c("Cage_Name","Genotype","Sex","Num_Mice","Num_Pups","Has_Pups",
                 "Status","Parent_Cage","Genotyping_Done","Genotype_Confirmed","Setup_Date"),
             drop = FALSE]

    # clickable parent cage link
    df$Parent_Cage <- ifelse(
      nzchar(df$Parent_Cage),
      sprintf("<a href='#' class='parent-cage' data-cage='%s'>%s</a>", df$Parent_Cage, df$Parent_Cage),
      ""
    )

    # add leading blank col for select checkbox
    df2 <- cbind(` ` = "", df)

    datatable(
      df2,
      extensions = c("Buttons","Select"),
      selection = "multiple",
      escape = FALSE, rownames = FALSE,
      options = list(
        pageLength = 10,
        dom = 'Bfrtip',
        buttons = c('copy','csv','excel','pdf','print'),
        order = list(list(1,'asc')), # sort by Cage_Name (now at col 1)
        columnDefs = list(
          list(orderable = FALSE, className = 'select-checkbox', targets = 0)
        ),
        select = list(style = 'multi', selector = 'td:first-child'),
        server = FALSE
      ),
      callback = JS("
        table.on('click', 'a.parent-cage', function(e){
          e.preventDefault();
          var cage = $(this).data('cage');
          Shiny.setInputValue('parent_cage_click', cage, {priority: 'event'});
        });
        // open cage when clicking its name (col 1)
        table.on('click', 'td:nth-child(2)', function(){
          var data = table.row(this).data();
          Shiny.setInputValue('cage_click', data[1], {priority:'event'});
        });
        // send selected Cage_Name keys to Shiny
        table.on('select.dt deselect.dt', function(e, dt){
          var cages = dt.rows({selected:true}).data().toArray().map(function(r){ return r[1]; });
          Shiny.setInputValue('cage_selected_keys_cv', cages, {priority:'event'});
        });
      ")
    ) %>%
      formatStyle('Cage_Name', cursor='pointer', color='#2196f3', fontWeight='bold') %>%
      formatStyle('Status',
        backgroundColor = styleEqual(
          c('Breeding','Holding','Experiment','Active','Retired'),
          c('#ffeb3b','#2196f3','#ff9800','#4caf50','#9e9e9e')
        )
      ) %>%
      formatStyle('Parent_Cage', cursor='pointer', color='#ff5722', fontWeight='bold') %>%
      formatStyle('Genotyping_Done',
        backgroundColor = styleEqual(c(TRUE,FALSE), c('#4caf50','#f44336'))
      ) %>%
      formatStyle('Genotype_Confirmed',
        backgroundColor = styleEqual(c(TRUE,FALSE), c('#4caf50','#f44336'))
      ) %>%
      formatStyle('Has_Pups',
        backgroundColor = styleEqual(c('Yes','No'), c('#e8f5e9','#ffebee')),
        color = styleEqual(c('Yes','No'), c('#1b5e20','#b71c1c')),
        fontWeight = 'bold'
      )
  })
  
  # Upload preview and issues outputs
  output$import_preview <- DT::renderDT({
    df <- uploaded_df()
    if (is.null(df)) return(DT::datatable(data.frame(Upload = "No file yet"), options = list(dom='t')))
    DT::datatable(df, extensions = "Buttons",
                  options = list(pageLength = 10, dom = 'Bfrtip',
                                 buttons = c('copy','csv','excel','pdf','print'),
                                 server = FALSE),
                  rownames = FALSE)
  })
  
  output$import_issues <- renderUI({
    iss <- upload_issues()
    if (length(iss$errors) == 0 && length(iss$warnings) == 0) return(HTML("<span style='color:#2e7d32'>No issues.</span>"))
    html <- ""
    if (length(iss$errors)) {
      html <- paste0(html, "<div style='color:#c62828'><b>Errors:</b><ul>",
                     paste(sprintf("<li>%s</li>", iss$errors), collapse = ""),
                     "</ul></div>")
    }
    if (length(iss$warnings)) {
      html <- paste0(html, "<div style='color:#ef6c00'><b>Warnings:</b><ul>",
                     paste(sprintf("<li>%s</li>", iss$warnings), collapse = ""),
                     "</ul></div>")
    }
    HTML(html)
  })
  
  # Gene table output
  output$gene_table <- DT::renderDT({
    df <- gene_catalog()
    df2 <- cbind(` ` = "", df)  # leading checkbox col
    datatable(
      df2,
      selection = "multiple",
      rownames = FALSE,
      extensions = c("Buttons","Select"),
      options = list(
        dom = 'Bfrtip',
        buttons = c('copy','csv','excel','pdf','print'),
        pageLength = 10,
        columnDefs = list(list(orderable = FALSE, className = 'select-checkbox', targets = 0)),
        select = list(style = 'single', selector = 'td:first-child')
      ),
      callback = JS("
        // Click the Gene name cell to load it into the editor
        table.on('click', 'td:nth-child(2)', function(){
          var data = table.row(this).data();
          if (data && data.length > 1) {
            Shiny.setInputValue('gene_row_clicked', data[1], {priority:'event'}); // Gene string
          }
        });

        // NEW: whenever selection changes, send selected Gene names
        table.on('select.dt deselect.dt', function(e, dt){
          var genes = dt.rows({selected:true}).data().toArray().map(function(r){ return r[1]; });
          Shiny.setInputValue('gene_selected_keys', genes, {priority:'event'});
        });
      ")
    ) %>%
      formatStyle('Gene', cursor='pointer', color='#2196f3', fontWeight='bold')
  })
  
  # Mouse filter UI outputs
  output$mouse_filter_gene_ui <- renderUI({
    genes <- unique(gene_catalog()$Gene)
    selectInput("mouse_filter_gene", "Gene", choices = c("Any", genes), selected = "Any")
  })
  
  output$mouse_filter_symbol_ui <- renderUI({
    g <- input$mouse_filter_gene %||% "Any"
    if (identical(g, "Any")) {
      selectInput("mouse_filter_symbol", "Genotype symbol", choices = "Any", selected = "Any")
    } else {
      syms <- get_gene_syms(g)
      selectInput("mouse_filter_symbol", "Genotype symbol", choices = c("Any", syms), selected = "Any")
    }
  })
  
  # --- Cage filters (Home) ---
  output$cage_filter_gene_ui <- renderUI({
    genes <- unique(gene_catalog()$Gene)
    selectInput("cage_filter_gene", "Gene", choices = c("Any", genes), selected = "Any")
  })
  output$cage_filter_symbol_ui <- renderUI({
    g <- input$cage_filter_gene %||% "Any"
    if (identical(g, "Any")) {
      selectInput("cage_filter_symbol", "Genotype symbol", choices = "Any", selected = "Any")
    } else {
      syms <- get_gene_syms(g)
      selectInput("cage_filter_symbol", "Genotype symbol", choices = c("Any", syms), selected = "Any")
    }
  })
  
  # --- Cage filters (Cage View list) ---
  output$cagecv_filter_gene_ui <- renderUI({
    genes <- unique(gene_catalog()$Gene)
    selectInput("cagecv_filter_gene", "Gene", choices = c("Any", genes), selected = "Any")
  })
  output$cagecv_filter_symbol_ui <- renderUI({
    g <- input$cagecv_filter_gene %||% "Any"
    if (identical(g, "Any")) {
      selectInput("cagecv_filter_symbol", "Genotype symbol", choices = "Any", selected = "Any")
    } else {
      syms <- get_gene_syms(g)
      selectInput("cagecv_filter_symbol", "Genotype symbol", choices = c("Any", syms), selected = "Any")
    }
  })
  
  # Keep cage choices fresh
  observe({
    cages <- sort(unique(cage_data_reactive()$Cage_Name))
    updateSelectizeInput(session, "mouse_filter_cage", choices = cages, server = TRUE)
  })
  
  # Reset button
  observeEvent(input$mouse_filters_reset, {
    updateTextInput(session, "mouse_search_text", value = "")
    updateSelectInput(session, "mouse_filter_gene", selected = "Any")
    updateSelectInput(session, "mouse_filter_symbol", selected = "Any")
    updateSelectInput(session, "mouse_filter_sex", selected = "Any")
    updateSelectInput(session, "mouse_filter_status", selected = "Any")
    updateSelectizeInput(session, "mouse_filter_cage", selected = character(0))
    updateDateRangeInput(session, "mouse_filter_dob", start = Sys.Date() - 90, end = Sys.Date())
    updateCheckboxInput(session, "hide_pups_mice", value = TRUE)
    updateCheckboxInput(session, "only_genotyped_mice", value = FALSE)
  })
  
  # Cage filter reset handlers
  observeEvent(input$cage_filters_reset, {
    updateSelectInput(session, "cage_filter_gene", selected = "Any")
    updateSelectInput(session, "cage_filter_symbol", selected = "Any")
    updateSelectInput(session, "cage_filter_sex", selected = "Any")
    updateSelectInput(session, "cage_filter_status", selected = "Any")
    updateSelectInput(session, "cage_filter_has_pups", selected = "Any")
    updateSelectInput(session, "cage_filter_genotyped", selected = "Any")
    updateSelectInput(session, "cage_filter_confirmed", selected = "Any")
    updateDateRangeInput(session, "cage_filter_setup", start = Sys.Date()-60, end = Sys.Date())
  })
  
  observeEvent(input$cagecv_filters_reset, {
    updateSelectInput(session, "cagecv_filter_gene", selected = "Any")
    updateSelectInput(session, "cagecv_filter_symbol", selected = "Any")
    updateSelectInput(session, "cagecv_filter_sex", selected = "Any")
    updateSelectInput(session, "cagecv_filter_status", selected = "Any")
    updateSelectInput(session, "cagecv_filter_has_pups", selected = "Any")
    updateSelectInput(session, "cagecv_filter_genotyped", selected = "Any")
    updateSelectInput(session, "cagecv_filter_confirmed", selected = "Any")
    updateDateRangeInput(session, "cagecv_filter_setup", start = Sys.Date()-60, end = Sys.Date())
  })
  
  # Filtered mice dataset
  filtered_mice <- reactive({
    df <- mice_data_reactive()

    # ensure columns exist
    need <- c("Mouse_ID","Cage","Sex","DOB","Genotype","Status","Notes","Mother_ID","Father_ID")
    for (nm in setdiff(need, names(df))) df[[nm]] <- NA

    # Free-text search
    q <- gsub("^\\s+|\\s+$", "", input$mouse_search_text %||% "")
    if (nzchar(q)) {
      keep <- grepl(q, df$Mouse_ID, TRUE) | grepl(q, df$Notes, TRUE) | grepl(q, df$Genotype, TRUE)
      df <- df[keep, , drop = FALSE]
    }

    # Gene / symbol filter
    g <- input$mouse_filter_gene %||% "Any"
    s <- input$mouse_filter_symbol %||% "Any"
    if (!identical(g, "Any")) {
      sym_here <- vapply(df$Genotype, extract_gene_symbol, character(1), gene = g)
      df <- df[!is.na(sym_here), , drop = FALSE]
      if (!identical(s, "Any")) df <- df[sym_here == s, , drop = FALSE]
    }

    # Sex
    sx <- input$mouse_filter_sex %||% "Any"
    if (!identical(sx, "Any")) df <- df[df$Sex == sx, , drop = FALSE]

    # Status
    st <- input$mouse_filter_status %||% "Any"
    if (!identical(st, "Any")) df <- df[tolower(df$Status) == tolower(st), , drop = FALSE]

    # Cage (multi)
    cs <- input$mouse_filter_cage %||% character(0)
    if (length(cs)) df <- df[df$Cage %in% cs, , drop = FALSE]

    # DOB range
    rng <- input$mouse_filter_dob
    if (!is.null(rng) && length(rng) == 2) {
      dts <- toDate_vec(df$DOB)
      df <- df[is.na(dts) | (dts >= as.Date(rng[1]) & dts <= as.Date(rng[2])), , drop = FALSE]
    }

    # Hide pups
    if (isTRUE(input$hide_pups_mice)) {
      df <- df[!(tolower(df$Status) == "pup"), , drop = FALSE]
    }

    # Only show rows that have at least one Gene(symbol) token, e.g., Cx3cr1(+/-)
    if (isTRUE(input$only_genotyped_mice)) {
      df <- df[grepl("\\w+\\([^)]*\\)", df$Genotype %||% ""), , drop = FALSE]
    }

    df
  })
  
  # Export filtered mice
  output$dl_filtered_mice <- downloadHandler(
    filename = function() sprintf("mLIMS_filtered_mice_%s.csv", Sys.Date()),
    content = function(file) utils::write.csv(filtered_mice(), file, row.names = FALSE, na = "")
  )
  
  # Export filtered cages (Home)
  output$dl_filtered_cages <- downloadHandler(
    filename = function() sprintf("mLIMS_filtered_cages_%s.csv", Sys.Date()),
    content  = function(file) {
      df <- filtered_cages_dashboard()
      df$Has_Pups <- ifelse(df$Num_Pups > 0, "Yes", "No")
      utils::write.csv(df, file, row.names = FALSE, na = "")
    }
  )
  
  # Export filtered cages (Cage View)
  output$dl_filtered_cages_cv <- downloadHandler(
    filename = function() sprintf("mLIMS_filtered_cages_cv_%s.csv", Sys.Date()),
    content  = function(file) {
      df <- filtered_cages_cv()
      df$Has_Pups <- ifelse(df$Num_Pups > 0, "Yes", "No")
      utils::write.csv(df, file, row.names = FALSE, na = "")
    }
  )
  
  # Home ▸ Cage Overview filtered data
  filtered_cages_dashboard <- reactive({
    df <- cage_data_reactive()

    # Search by name
    q <- gsub("^\\s+|\\s+$", "", input$cage_search %||% "")
    if (nzchar(q)) df <- df[grepl(q, df$Cage_Name, ignore.case = TRUE), , drop = FALSE]

    # Quick toggle: pending genotyping
    if (isTRUE(input$filter_genotyping)) df <- df[!df$Genotyping_Done, , drop = FALSE]

    # Compute Has_Pups for filtering
    df$Has_Pups <- ifelse(df$Num_Pups > 0, "Yes", "No")

    # Gene / symbol filter (uses your extract_gene_symbol)
    g <- input$cage_filter_gene   %||% "Any"
    s <- input$cage_filter_symbol %||% "Any"
    if (!identical(g, "Any")) {
      sym_here <- vapply(df$Genotype, extract_gene_symbol, character(1), gene = g)
      keep <- !is.na(sym_here)
      if (!identical(s, "Any")) keep <- keep & sym_here == s
      df <- df[keep, , drop = FALSE]
    }

    # Sex, Status
    sx <- input$cage_filter_sex    %||% "Any"
    st <- input$cage_filter_status %||% "Any"
    if (!identical(sx, "Any")) df <- df[df$Sex == sx, , drop = FALSE]
    if (!identical(st, "Any")) df <- df[df$Status == st, , drop = FALSE]

    # Has pups
    hp <- input$cage_filter_has_pups %||% "Any"
    if (hp %in% c("Yes","No")) df <- df[df$Has_Pups == hp, , drop = FALSE]

    # Genotyping flags
    gd <- input$cage_filter_genotyped %||% "Any"
    gc <- input$cage_filter_confirmed %||% "Any"
    if (gd %in% c("Yes","No")) df <- df[df$Genotyping_Done == (gd == "Yes"), , drop = FALSE]
    if (gc %in% c("Yes","No")) df <- df[df$Genotype_Confirmed == (gc == "Yes"), , drop = FALSE]

    # Setup date range
    rng <- input$cage_filter_setup
    if (!is.null(rng) && length(rng) == 2) {
      dts <- toDate_vec(df$Setup_Date)
      df <- df[is.na(dts) | (dts >= as.Date(rng[1]) & dts <= as.Date(rng[2])), , drop = FALSE]
    }

    df
  })
  
  # 🔎 Cage View (list) filtered data
  filtered_cages_cv <- reactive({
    df <- cage_data_reactive()

    # Search
    q <- gsub("^\\s+|\\s+$", "", input$cage_search_cv %||% "")
    if (nzchar(q)) df <- df[grepl(q, df$Cage_Name, ignore.case = TRUE), , drop = FALSE]

    # Quick toggle
    if (isTRUE(input$filter_genotyping_cv)) df <- df[!df$Genotyping_Done, , drop = FALSE]

    # Has_Pups for filtering
    df$Has_Pups <- ifelse(df$Num_Pups > 0, "Yes", "No")

    # Gene/symbol
    g <- input$cagecv_filter_gene   %||% "Any"
    s <- input$cagecv_filter_symbol %||% "Any"
    if (!identical(g, "Any")) {
      sym_here <- vapply(df$Genotype, extract_gene_symbol, character(1), gene = g)
      keep <- !is.na(sym_here)
      if (!identical(s, "Any")) keep <- keep & sym_here == s
      df <- df[keep, , drop = FALSE]
    }

    # Sex, Status
    sx <- input$cagecv_filter_sex    %||% "Any"
    st <- input$cagecv_filter_status %||% "Any"
    if (!identical(sx, "Any")) df <- df[df$Sex == sx, , drop = FALSE]
    if (!identical(st, "Any")) df <- df[df$Status == st, , drop = FALSE]

    # Has pups
    hp <- input$cagecv_filter_has_pups %||% "Any"
    if (hp %in% c("Yes","No")) df <- df[df$Has_Pups == hp, , drop = FALSE]

    # Genotyping flags
    gd <- input$cagecv_filter_genotyped %||% "Any"
    gc <- input$cagecv_filter_confirmed %||% "Any"
    if (gd %in% c("Yes","No")) df <- df[df$Genotyping_Done == (gd == "Yes"), , drop = FALSE]
    if (gc %in% c("Yes","No")) df <- df[df$Genotype_Confirmed == (gc == "Yes"), , drop = FALSE]

    # Setup date
    rng <- input$cagecv_filter_setup
    if (!is.null(rng) && length(rng) == 2) {
      dts <- toDate_vec(df$Setup_Date)
      df <- df[is.na(dts) | (dts >= as.Date(rng[1]) & dts <= as.Date(rng[2])), , drop = FALSE]
    }

    df
  })
  
  # Upload modal handlers
  observeEvent(input$shortcut_upload, {
    showModal(uploadModalUI(!is.null(mlims_cols())))
  })
  
  observeEvent(input$open_upload, {
    showModal(uploadModalUI(!is.null(mlims_cols())))
  })
  
  # Optional: keep keyboard shortcuts to trigger the download buttons
  observeEvent(input$shortcut_download_all, {
    runjs("document.getElementById('dl_all').click();")
  })
  
  observeEvent(input$shortcut_export_cage, {
    if (!is.null(current_cage())) {
      runjs("document.getElementById('dl_cage').click();")
    }
  })
  
  observeEvent(input$shortcut_print_all, {
    showNotification("🖨️ Opening print preview for all cages...", type = "message")
  })
  
  # From inside an opened cage (button and keyboard shortcut both land here)
  observeEvent(input$add_mouse_to_cage, {
    req(current_cage())
    showAddMouseModal(default_cage = current_cage())
  })
  
  observeEvent(input$shortcut_add_mouse, {
    if (!is.null(current_cage())) {
      showAddMouseModal(default_cage = current_cage())
    } else {
      showAddMouseModal(default_cage = NULL)  # global add (choose/create cage)
    }
  })
  
  # From the Mouse Records page (global add)
  observeEvent(input$add_mouse, {
    showAddMouseModal(default_cage = NULL)  # user picks or creates a cage here
  })
  
  # One save handler that creates/updates the target cage automatically
  observeEvent(input$save_new_mouse_unified, {
    cage_name <- gsub("^\\s+|\\s+$", "", input$new_mouse_cage_select %||% "")
    req(nzchar(cage_name), input$new_mouse_id_unified, input$new_mouse_sex_unified)

    all <- mice_data_reactive()

    # Ensure union of columns exists
    new_cols <- c("Mouse_ID","Cage","Sex","DOB","Genotype","Status","Notes",
                  "Weight","Health_Status","Strain","Source","Lab_ID","Experiment_Group",
                  "Mother_ID","Father_ID")
    for (nm in setdiff(new_cols, names(all))) all[[nm]] <- NA

    row <- data.frame(
      Mouse_ID         = input$new_mouse_id_unified,
      Cage             = cage_name,
      Sex              = input$new_mouse_sex_unified,
      DOB              = input$new_mouse_dob_unified,
      Genotype         = input$new_mouse_genotype_unified %||% "Unknown",
      Status           = input$new_mouse_status_unified %||% "Active",
      Notes            = input$new_mouse_notes_unified %||% "",
      Weight           = input$new_mouse_weight_unified,
      Health_Status    = input$new_mouse_health_unified %||% "",
      Strain           = input$new_mouse_strain_unified %||% "",
      Source           = input$new_mouse_source_unified %||% "",
      Lab_ID           = input$new_mouse_labid_unified %||% "",
      Experiment_Group = input$new_mouse_group_unified %||% "",
      Mother_ID        = NA_character_,
      Father_ID        = NA_character_,
      stringsAsFactors = FALSE
    )
    row <- row[, names(all), drop = FALSE]

    # Append mouse
    all <- rbind(all, row)
    mice_data_reactive(all)

    # Ensure cage exists (or minimally update it). If it's a new name, this creates the cage.
    add_or_update_cage(
      cage_name = cage_name,
      genotype  = "Unknown",
      sex       = "Mixed",     # will correct below based on mice present
      status    = "Active",
      parent    = "",
      setup     = format(Sys.Date(), "%m/%d")
    )

    recalc_cage_counts()
    set_cage_sex_from_mice(c(cage_name))  # auto-set cage sex from its mice

    # Log the action
    log_event(
      "Add Mouse",
      cage = cage_name,
      mouse = input$new_mouse_id_unified,
      details = paste("Sex:", input$new_mouse_sex_unified,
                      "| DOB:", input$new_mouse_dob_unified %||% "",
                      "| Genotype:", input$new_mouse_genotype_unified %||% "Unknown",
                      "| Status:", input$new_mouse_status_unified %||% "Active")
    )

    removeModal()
    showNotification(paste("✅ Mouse added to cage", cage_name), type = "message")
  })
  

  
  observeEvent(input$shortcut_print_cage, {
    if (!is.null(current_cage())) {
      showNotification(paste("🖨️ Opening print preview for cage", current_cage(), "..."), type = "message")
    }
  })
  
  observeEvent(input$shortcut_edit_mouse, {
    if (!is.null(current_mouse())) {
      showModal(modalDialog(
        title = paste("Edit Mouse:", current_mouse()),
        div(
          textInput("edit_mouse_id", "Mouse ID:", value = current_mouse()),
          selectInput("edit_mouse_sex", "Sex:", choices = c("M", "F")),
          textInput("edit_mouse_dob", "Date of Birth:"),
          selectInput("edit_mouse_genotype", "Genotype:", choices = c("+/+", "+/-", "-/-")),
          textInput("edit_mouse_notes", "Notes:")
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("save_edited_mouse", "Save", class = "action-btn btn-view")
        )
      ))
    }
  })
  
  observeEvent(input$save_edited_mouse, {
    req(current_mouse())
    all <- mice_data_reactive()
    i <- which(all$Mouse_ID == current_mouse())[1]
    if (!is.na(i)) {
      all$Mouse_ID[i] <- input$edit_mouse_id %||% all$Mouse_ID[i]
      all$Sex[i]      <- input$edit_mouse_sex %||% all$Sex[i]
      all$DOB[i]      <- input$edit_mouse_dob %||% all$DOB[i]
      all$Genotype[i] <- input$edit_mouse_genotype %||% all$Genotype[i]
      all$Notes[i]    <- input$edit_mouse_notes %||% all$Notes[i]
      mice_data_reactive(all)
      recalc_cage_counts()
      set_cage_sex_from_mice(all$Cage[i])
      removeModal()
      showNotification("✅ Mouse updated.", type="message")
    }
  })
  
  observeEvent(input$shortcut_back_home, {
    # If we're in Cage View with a selected cage, go back to cage list
    if (identical(input$sidebar, "cage_view") && !is.null(current_cage())) {
      current_cage(NULL)
    } else {
      # Otherwise go back to home
      updateTabItems(session, "sidebar", "dashboard")
      current_cage(NULL)
      current_mouse(NULL)
    }
  })
  
  # Handle parent cage clicks
  observeEvent(input$parent_cage_click, {
    req(input$parent_cage_click)
    current_cage(input$parent_cage_click)
    current_mouse(NULL)
    updateTabItems(session, "sidebar", "cage_view")
  })
  
  # Quick action handlers
  
  observeEvent(input$add_pup_to_cage, {
    req(current_cage())

    all <- mice_data_reactive()
    in_cage <- all[all$Cage == current_cage(), , drop = FALSE]
    # adults (exclude pups)
    adults  <- in_cage[!is_pup_row(in_cage), , drop = FALSE]
    moms    <- adults$Mouse_ID[adults$Sex == "F"]
    dads    <- adults$Mouse_ID[adults$Sex == "M"]

    showModal(modalDialog(
      title = paste("Add Pup to", current_cage()),
      div(
        numericInput("addpup_n", "Number of pups:", value = 2, min = 1, max = 30),
        selectInput("addpup_sex", "Sex distribution:", c("Mixed","M","F"), selected="Mixed"),
        textInput("addpup_prefix", "Pup ID prefix:",
                  value = paste0(gsub("\\W","", current_cage()), "-P")),
        dateInput("addpup_dob", "DOB:", value = Sys.Date()),
        textInput("addpup_gt", "Genotype:", value = "Unknown"),

        tags$hr(),
        h4("Parents (optional)"),
        # Mother picker
        if (length(moms)) {
          selectizeInput(
            "addpup_mother_sel", "Mother (♀) in this cage:",
            choices  = moms,
            selected = if (length(moms) == 1) moms[1] else NULL,
            options  = list(placeholder = "Select a female or leave blank")
          )
        } else {
          tags$small(style="color:#777;", "No adult females available in this cage.")
        },
        # Father picker
        if (length(dads)) {
          selectizeInput(
            "addpup_father_sel", "Father (♂) in this cage:",
            choices  = dads,
            selected = if (length(dads) == 1) dads[1] else NULL,
            options  = list(placeholder = "Select a male or leave blank")
          )
        } else {
          tags$small(style="color:#777;", "No adult males available in this cage.")
        },

        tags$hr(),
        textInput("addpup_notes", "Notes:", value = "Born")
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_addpup", "Add Pup(s)", class="action-btn btn-view")
      )
    ))
  })
  
  observeEvent(input$save_addpup, {
    req(current_cage(), input$addpup_n, input$addpup_prefix)
    n <- as.integer(input$addpup_n); if (is.na(n) || n < 1) return()

    all <- mice_data_reactive()
    # ensure parent columns exist
    if (!"Mother_ID" %in% names(all)) all$Mother_ID <- NA_character_
    if (!"Father_ID" %in% names(all)) all$Father_ID <- NA_character_

    # Recompute valid parent pools for safety
    in_cage <- all[all$Cage == current_cage(), , drop = FALSE]
    adults  <- in_cage[!is_pup_row(in_cage), , drop = FALSE]
    moms_ok <- adults$Mouse_ID[adults$Sex == "F"]
    dads_ok <- adults$Mouse_ID[adults$Sex == "M"]

    mom <- input$addpup_mother_sel %||% ""
    dad <- input$addpup_father_sel %||% ""

    # Validate selections if provided
    if (nzchar(mom) && !(mom %in% moms_ok)) {
      showNotification("Selected mother is not an adult female in this cage.", type="error"); return()
    }
    if (nzchar(dad) && !(dad %in% dads_ok)) {
      showNotification("Selected father is not an adult male in this cage.", type="error"); return()
    }

    ids <- make_unique_ids(input$addpup_prefix, n, all$Mouse_ID)
    sex <- switch(input$addpup_sex,
                  "M"="M","F"="F",
                  Mixed = ifelse((seq_len(n) %% 2)==1, "M", "F"))

    new_rows <- data.frame(
      Mouse_ID  = ids,
      Cage      = current_cage(),
      Sex       = sex,
      DOB       = format(input$addpup_dob, "%Y-%m-%d"),
      Genotype  = input$addpup_gt %||% "Unknown",
      Status    = "Pup",
      Notes     = input$addpup_notes %||% "Born",
      Mother_ID = if (nzchar(mom)) mom else NA_character_,
      Father_ID = if (nzchar(dad)) dad else NA_character_,
      stringsAsFactors = FALSE
    )

    mice_data_reactive(rbind(all, new_rows))
    recalc_cage_counts()
    
    # Log the action
    log_event(
      "Add Pup(s)",
      cage = current_cage(),
      count = n,
      details = paste(
        "DOB:", format(input$addpup_dob, "%Y-%m-%d"),
        "| Sex:", input$addpup_sex,
        "| Prefix:", input$addpup_prefix,
        if (nzchar(input$addpup_mother_sel %||% "")) paste("| Mother:", input$addpup_mother_sel) else "",
        if (nzchar(input$addpup_father_sel %||% "")) paste("| Father:", input$addpup_father_sel) else ""
      )
    )
    
    removeModal()
    showNotification(
      paste0("✅ Added ", n, " pup(s) to ", current_cage(),
             if (nzchar(mom) || nzchar(dad)) {
               paste0(" (Mother: ", ifelse(nzchar(mom), mom, "—"),
                      ", Father: ", ifelse(nzchar(dad), dad, "—"), ")")
             } else ""), 
      type="message"
    )
  })
  
  observeEvent(input$mark_weaned, {
    req(current_cage())
    n_in_cage <- nrow(mice_data_reactive()[mice_data_reactive()$Cage == current_cage(), ])

    showModal(modalDialog(
      title = paste("Wean / Split from", current_cage()),
      size = "l",
      easyClose = TRUE,
      div(
        p("Optionally select pups in the table first, OR specify counts below."),
        dateInput("weaning_date", "Weaning Date:", value = Sys.Date()),

        # Use selection (if any)
        checkboxInput("wean_use_selected",
                      "Use selected pups from table (if any)",
                      value = n_in_cage > 0),

        # Global ID prefix (used for any created pups)
        textInput("pup_prefix", "Pup ID prefix:",
                  value = paste0(gsub("\\W","", current_cage()), "-P")),

        checkboxInput("wean_fill_shortfall",
                      "If not enough pups exist, create new pups to fill requested numbers",
                      value = TRUE),
        selectInput("post_wean_status",
                    "Status for weaned mice:",
                    c("Juvenile","Active","Holding"),
                    selected = "Juvenile"),

        tags$hr(),
        numericInput("dest_count", "Number of destination cages:", value = 2, min = 1, max = 10),
        uiOutput("split_destinations"),

        textInput("weaning_notes", "Notes:", value = "Normal weaning")
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_weaning_multi", "Split / Wean", class = "action-btn btn-edit")
      )
    ))
  })
  

  
  output$split_destinations <- renderUI({
    n <- input$dest_count %||% 1
    tagList(lapply(seq_len(n), function(i) {
      fluidRow(
        column(4, textInput(paste0("dest_name_", i),
                            paste0("Cage ", i, " name:"),
                            value = paste0(current_cage(), "-W", i))),
        column(2, numericInput(paste0("dest_pup_n_", i),
                               "Pup number", value = 2, min = 0, max = 50)),
        column(2, selectInput(paste0("dest_pup_sex_", i),
                              "Pup sex", c("Mixed","M","F"), selected = "Mixed")),
        column(4, textInput(paste0("dest_genotype_", i),
                            "Genotype", value = "Unknown"))
      )
    }))
  })
  
  observeEvent(input$save_weaning_multi, {
    req(current_cage())

    # --- collect destination cages + requested counts/sex ---
    n <- input$dest_count %||% 1
    dests <- lapply(seq_len(n), function(i) {
      list(
        name = input[[paste0("dest_name_", i)]],
        pup_n = as.integer(input[[paste0("dest_pup_n_", i)]]) %||% 0,
        pup_sex_mode = input[[paste0("dest_pup_sex_", i)]] %||% "Mixed",
        gt = input[[paste0("dest_genotype_", i)]] %||% "Unknown"
      )
    })
    if (any(vapply(dests, function(d) !nzchar(d$name), logical(1)))) {
      showNotification("All destination cages need a name.", type="error"); return()
    }
    total_req <- sum(vapply(dests, `[[`, integer(1), "pup_n"))
    if (total_req < 1) {
      showNotification("Set 'Pup number' > 0 for at least one destination.", type="error"); return()
    }

    # Ensure/prepare destination cages with minimal metadata
    lapply(dests, function(d) add_or_update_cage(
      d$name, genotype = d$gt, sex = "Mixed",
      status = "Holding", parent = current_cage(),
      setup = format(input$weaning_date, "%m/%d")
    ))

    all_mice <- mice_data_reactive()
    in_src   <- all_mice[all_mice$Cage == current_cage(), , drop=FALSE]

    # pool of existing pups in the source cage
    pool <- in_src[is_pup_row(in_src), , drop = FALSE]

    # Which rows did the user select (if any)?
    sel_idx <- if (isTRUE(input$use_mlims_template))
                 input$mice_in_cage_table_mlims_rows_selected
               else
                 input$mice_in_cage_table_rows_selected
    wants_selection <- isTRUE(input$wean_use_selected)
    selected <- if (wants_selection && length(sel_idx)) in_src[sel_idx, , drop=FALSE] else in_src[0, , drop=FALSE]
    # keep only pups from the explicit selection
    if (nrow(selected)) selected <- selected[is_pup_row(selected), , drop = FALSE]

    # destination capacities and counters
    caps <- vapply(dests, `[[`, integer(1), "pup_n")
    assigned_counts <- integer(length(caps))
    names(caps) <- names(assigned_counts) <- vapply(dests, `[[`, "", "name")

    # helper: take up to k from df by preferred sex, falling back to any
    take_from_pool <- function(df, k, sex_pref) {
      if (k <= 0 || !nrow(df)) return(list(taken=df[0,], rest=df))
      if (sex_pref %in% c("M","F")) {
        want <- df[df$Sex == sex_pref, , drop=FALSE]
        if (nrow(want) >= k) {
          taken <- head(want, k)
          rest  <- df[!(df$Mouse_ID %in% taken$Mouse_ID), , drop=FALSE]
          return(list(taken=taken, rest=rest))
        } else {
          # take what we can, then top up with any sex
          taken1 <- want
          still  <- k - nrow(taken1)
          any_left <- df[!(df$Mouse_ID %in% taken1$Mouse_ID), , drop=FALSE]
          taken2 <- head(any_left, still)
          taken <- rbind(taken1, taken2)
          rest  <- any_left[!(any_left$Mouse_ID %in% taken2$Mouse_ID), , drop=FALSE]
          return(list(taken=taken, rest=rest))
        }
      } else {
        # Mixed: alternate if possible
        mm <- df[df$Sex == "M", , drop=FALSE]
        ff <- df[df$Sex == "F", , drop=FALSE]
        picked <- df[0, , drop=FALSE]
        turn <- "M"
        while (nrow(picked) < k && (nrow(mm) > 0 || nrow(ff) > 0)) {
          if (turn == "M" && nrow(mm) > 0) { picked <- rbind(picked, mm[1,]); mm <- mm[-1,,drop=FALSE]; turn <- "F"; next }
          if (turn == "F" && nrow(ff) > 0) { picked <- rbind(picked, ff[1,]); ff <- ff[-1,,drop=FALSE]; turn <- "M"; next }
          # fallback: whichever remains
          if (nrow(mm) > 0) { picked <- rbind(picked, mm[1,]); mm <- mm[-1,,drop=FALSE]; next }
          if (nrow(ff) > 0) { picked <- rbind(picked, ff[1,]); ff <- ff[-1,,drop=FALSE]; next }
        }
        rest <- rbind(mm, ff)
        return(list(taken=picked, rest=rest))
      }
    }

    moved <- in_src[0, , drop=FALSE]

    # 1) If user explicitly selected pups, use them first (respect cage caps in round-robin)
    if (nrow(selected)) {
      d_names <- names(caps)
      j <- 1
      for (k in seq_len(nrow(selected))) {
        # find next cage with remaining capacity
        tried <- 0
        while (caps[j] - assigned_counts[j] <= 0 && tried < length(caps)) {
          j <- if (j == length(caps)) 1 else j + 1
          tried <- tried + 1
        }
        if (caps[j] - assigned_counts[j] > 0) {
          moved <- rbind(moved, selected[k, , drop=FALSE])
          moved$Cage[nrow(moved)] <- d_names[j]
          assigned_counts[j] <- assigned_counts[j] + 1
        }
        j <- if (j == length(caps)) 1 else j + 1
      }
      # remove those from pool too so they aren't double-counted
      if (nrow(moved)) {
        pool <- pool[!(pool$Mouse_ID %in% moved$Mouse_ID), , drop=FALSE]
      }
      # after you've built `moved` from selected rows
      if (nrow(moved)) moved$Status <- input$post_wean_status
    }

    # 2) Auto-pick remaining required pups from the pool to satisfy each destination
    for (i in seq_along(dests)) {
      need <- caps[i] - assigned_counts[i]
      if (need <= 0) next
      pick <- take_from_pool(pool, need, dests[[i]]$pup_sex_mode)
      got  <- pick$taken; pool <- pick$rest
      if (nrow(got)) {
        got$Cage   <- dests[[i]]$name
        got$Status <- input$post_wean_status
        moved <- rbind(moved, got)
        assigned_counts[i] <- assigned_counts[i] + nrow(got)
      }
    }

    # 3) Optional: create new pups to fill shortfall (if user allowed)
    created_total <- 0
    if (isTRUE(input$wean_fill_shortfall)) {
      make_ids <- function(n) make_unique_ids(input$pup_prefix, n, all_mice$Mouse_ID)
      for (i in seq_along(dests)) {
        need <- caps[i] - assigned_counts[i]
        if (need <= 0) next
        d <- dests[[i]]
        new_ids <- make_ids(need)
        new_sex <- if (d$pup_sex_mode == "M") rep("M", need) else if (d$pup_sex_mode == "F") rep("F", need)
                   else if (need>0) ifelse((seq_len(need) %% 2)==1, "M", "F")
        new_rows <- data.frame(
          Mouse_ID = new_ids,
          Cage     = d$name,
          Sex      = new_sex,
          DOB      = format(input$weaning_date, "%Y-%m-%d"),
          Genotype = d$gt,
          Status   = input$post_wean_status,
          Notes    = paste("Weaned", input$weaning_notes %||% ""),
          stringsAsFactors = FALSE
        )
        all_mice <- rbind(all_mice, new_rows)
        created_total <- created_total + need
      }
    } else {
      # If not creating new, warn if short
      short <- caps - assigned_counts
      if (any(short > 0)) {
        showNotification(
          paste("Not enough pups available; unfilled:", 
                paste(names(short[short>0]), short[short>0], collapse = "; ")),
          type="warning"
        )
      }
    }

    # 4) Apply moves: remove moved-from-source rows then add back with new cage
    if (nrow(moved)) {
      all_mice <- all_mice[!(all_mice$Mouse_ID %in% moved$Mouse_ID), ]
      all_mice <- rbind(all_mice, moved)
    }

    mice_data_reactive(all_mice)
    recalc_cage_counts()
    set_cage_sex_from_mice(vapply(dests, `[[`, "", "name"))

    removeModal()
    showNotification(
      paste0("✅ Weaned/split into ", length(caps), " cage(s). ",
             nrow(moved), " moved, ",
             created_total, " created."),
      type="message"
    )

    current_cage(NULL)  # back to list so you can see new cages
  })
  
  observeEvent(input$split_cage, {
    req(current_cage())
    # how many in source (for hint only)
    n_in_cage <- nrow(mice_data_reactive()[mice_data_reactive()$Cage == current_cage(), ])
    showModal(modalDialog(
      title = paste("Split Cage:", current_cage()),
      size = "l",
      easyClose = TRUE,
      div(
        p("Option 1) Select rows in the table first, then open this dialog."),
        p("Option 2) Or specify counts below (will auto-pick from this cage)."),
        checkboxInput("split_use_selected", "Use selected rows if any", value = n_in_cage > 0),
        tags$hr(),
        numericInput("split_dest_count", "Number of destination cages:", value = 2, min = 1, max = 10),
        uiOutput("split_destinations_ui"),
        textInput("split_notes", "Notes:", value = "Cage split")
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_split_cage", "Split", class="action-btn btn-edit")
      )
    ))
  })
  
  output$split_destinations_ui <- renderUI({
    req(current_cage())
    n <- input$split_dest_count %||% 1
    tagList(lapply(seq_len(n), function(i) {
      fluidRow(
        column(6, textInput(paste0("split_dest_name_", i),
                            paste0("Cage ", i, " name:"),
                            value = paste0(current_cage(), "-S", i))),
        column(3, numericInput(paste0("split_dest_n_", i),
                               "Number of mice", value = 1, min = 0, max = 99)),
        column(3, selectInput(paste0("split_dest_sex_", i),
                              "Prefer sex", c("Any","M","F"), selected="Any"))
      )
    }))
  })
  
  observeEvent(input$save_split_cage, {
    req(current_cage())
    # gather destinations
    n <- input$split_dest_count %||% 1
    dests <- lapply(seq_len(n), function(i) {
      list(
        name = input[[paste0("split_dest_name_", i)]],
        want = as.integer(input[[paste0("split_dest_n_", i)]]) %||% 0,
        sex  = input[[paste0("split_dest_sex_", i)]] %||% "Any"
      )
    })
    if (any(vapply(dests, function(d) !nzchar(d$name), logical(1)))) {
      showNotification("All destination cages need a name.", type="error"); return()
    }
    total_req <- sum(vapply(dests, `[[`, integer(1), "want"))
    if (total_req < 1) { showNotification("Set 'Number of mice' > 0.", type="error"); return() }

    # ensure cages exist
    lapply(dests, function(d) add_or_update_cage(
      d$name, genotype = cage_data_reactive()$Genotype[cage_data_reactive()$Cage_Name==current_cage()][1] %||% "Unknown",
      sex = "Mixed", status = "Holding", parent = current_cage(),
      setup = format(Sys.Date(), "%m/%d")
    ))

    all <- mice_data_reactive()
    src <- all[all$Cage == current_cage(), , drop=FALSE]

    # selected rows?
    sel_idx <- if (isTRUE(input$use_mlims_template))
                 input$mice_in_cage_table_mlims_rows_selected
               else
                 input$mice_in_cage_table_rows_selected
    wants_sel <- isTRUE(input$split_use_selected)
    selected <- if (wants_sel && length(sel_idx)) src[sel_idx, , drop=FALSE] else src[0, , drop=FALSE]

    # assign selected first (round-robin into cages with remaining capacity)
    caps <- vapply(dests, `[[`, integer(1), "want")
    moved <- src[0, , drop=FALSE]
    if (nrow(selected)) {
      j <- 1; names_caps <- vapply(dests, `[[`, "", "name"); used <- integer(length(caps))
      for (k in seq_len(nrow(selected))) {
        tried <- 0
        while (caps[j] - used[j] <= 0 && tried < length(caps)) { j <- if (j==length(caps)) 1 else j+1; tried <- tried+1 }
        if (caps[j] - used[j] > 0) {
          moved <- rbind(moved, selected[k, , drop=FALSE])
          moved$Cage[nrow(moved)] <- names_caps[j]
          used[j] <- used[j] + 1
        }
        j <- if (j==length(caps)) 1 else j+1
      }
      # remove moved ones from pool
      src <- src[!(src$Mouse_ID %in% moved$Mouse_ID), , drop=FALSE]
      caps <- pmax(0L, caps - table(factor(moved$Cage, levels=names_caps)))
      caps <- as.integer(caps)
    }

    # auto-pick to fill remaining requests
    pick_pref <- function(df, k, sex){
      if (k<=0 || !nrow(df)) return(list(take=df[0,], rest=df))
      if (sex %in% c("M","F")) {
        pref <- df[df$Sex==sex, , drop=FALSE]
        if (nrow(pref)>=k) return(list(take=head(pref,k), rest=df[!(df$Mouse_ID %in% head(pref,k)$Mouse_ID),,drop=FALSE]))
        take1 <- pref; still <- k-nrow(pref)
        remain <- df[!(df$Mouse_ID %in% pref$Mouse_ID), , drop=FALSE]
        take2 <- head(remain, still)
        take <- rbind(take1, take2)
        rest <- remain[!(remain$Mouse_ID %in% take2$Mouse_ID), , drop=FALSE]
        return(list(take=take, rest=rest))
      } else {
        take <- head(df, k); rest <- df[-seq_len(nrow(take)), , drop=FALSE]; return(list(take=take, rest=rest))
      }
    }

    for (i in seq_along(dests)) {
      need <- caps[i]; if (need<=0) next
      pp <- pick_pref(src, need, dests[[i]]$sex)
      tk <- pp$take; src <- pp$rest
      if (nrow(tk)) { tk$Cage <- dests[[i]]$name; moved <- rbind(moved, tk) }
    }

    if (!nrow(moved)) { showNotification("No mice moved.", type="warning"); return() }

    # apply move
    all <- all[!(all$Mouse_ID %in% moved$Mouse_ID), ]
    all <- rbind(all, moved)
    mice_data_reactive(all)
    recalc_cage_counts()
    set_cage_sex_from_mice(vapply(dests, `[[`, "", "name"))

    removeModal()
    showNotification(paste0("✅ Split complete: moved ", nrow(moved), " mouse/mice into ",
                            length(dests), " cage(s)."), type="message")
  })
  
  observeEvent(input$edit_cage, {
    req(current_cage())
    cg <- cage_data_reactive()
    row <- cg[cg$Cage_Name == current_cage(), , drop=FALSE]
    showModal(modalDialog(
      title = paste("Edit Cage:", current_cage()),
      div(
        textInput("edit_cage_name", "Cage Name:", value = row$Cage_Name[1]),
        selectInput("edit_cage_status", "Status:",
                    c("Breeding","Holding","Experiment","Active","Retired"),
                    selected = row$Status[1]),
        textInput("edit_cage_genotype", "Genotype:", value = row$Genotype[1]),
        selectInput("edit_cage_sex", "Sex:", c("M","F","Mixed"), selected = row$Sex[1])
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_edit_cage", "Save", class = "action-btn btn-view")
      )
    ))
  })

  observeEvent(input$save_edit_cage, {
    cg <- cage_data_reactive()
    i <- which(cg$Cage_Name == current_cage())[1]
    if (!is.na(i)) {
      cg$Cage_Name[i] <- input$edit_cage_name
      cg$Status[i]    <- input$edit_cage_status
      cg$Genotype[i]  <- input$edit_cage_genotype
      cg$Sex[i]       <- input$edit_cage_sex
      cage_data_reactive(cg)
      current_cage(input$edit_cage_name)  # keep context if name changed
    }
    removeModal()
    showNotification("✅ Cage updated.", type="message")
  })

  # Gene catalog handlers
  observeEvent(input$add_gene_btn, {
    nm <- gsub("^\\s+|\\s+$", "", input$gene_name %||% "")
    sy <- normalize_syms(input$gene_symbols)
    if (!nzchar(nm) || !nzchar(sy)) {
      showNotification("Enter a gene name and at least one symbol.", type="warning"); return()
    }
    gc <- gene_catalog()
    if (nm %in% gc$Gene) {
      gc$Symbols[gc$Gene == nm] <- sy
    } else {
      gc <- rbind(gc, data.frame(Gene = nm, Symbols = sy, stringsAsFactors = FALSE))
    }
    gene_catalog(gc)
    last_gene(nm)
    showNotification(paste0("✅ Gene catalog updated: ", nm), type="message")
  })
  
  observeEvent(input$delete_gene_btn, {
    keys <- input$gene_selected_keys %||% character(0)
    gc <- gene_catalog()
    if (!length(keys)) { 
      showNotification("Select one or more genes to delete.", type="warning"); return() 
    }
    gene_catalog(gc[!(gc$Gene %in% keys), , drop = FALSE])
    if (!is.null(last_gene()) && last_gene() %in% keys) last_gene(NULL)
    showNotification(paste0("🗑️ Deleted: ", paste(keys, collapse = ", ")), type="message")
  })
  
  # Fill editor when a row is selected via checkbox
  observeEvent(input$gene_table_rows_selected, {
    sel <- input$gene_table_rows_selected
    if (length(sel) == 1) {
      gc <- gene_catalog()
      updateTextInput(session, "gene_name", value = gc$Gene[sel])
      updateTextInput(session, "gene_symbols", value = gc$Symbols[sel])
    }
  })
  
  # Fill editor when genes are selected via stable keys
  observeEvent(input$gene_selected_keys, {
    ks <- input$gene_selected_keys
    if (length(ks) == 1) {
      row <- subset(gene_catalog(), Gene == ks[1])
      updateTextInput(session, "gene_name", value = row$Gene[1])
      updateTextInput(session, "gene_symbols", value = row$Symbols[1])
    }
  })
  
  # Fill editor when you click the Gene name cell
  observeEvent(input$gene_row_clicked, {
    g <- input$gene_row_clicked
    gc <- gene_catalog()
    if (nrow(gc) && g %in% gc$Gene) {
      row <- gc[gc$Gene == g, , drop = FALSE]
      updateTextInput(session, "gene_name", value = row$Gene[1])
      updateTextInput(session, "gene_symbols", value = row$Symbols[1])
    }
  })
  
  # Handler to save a quick gene and then reopen Genotype modal
  observeEvent(input$quick_add_gene_save, {
    nm <- gsub("^\\s+|\\s+$", "", input$quick_gene_name %||% "")
    sy <- normalize_syms(input$quick_gene_syms)
    if (!nzchar(nm) || !nzchar(sy)) {
      showNotification("Enter a gene name and at least one symbol.", type="warning"); return()
    }
    gc <- gene_catalog()
    if (nm %in% gc$Gene) {
      gc$Symbols[gc$Gene == nm] <- sy
    } else {
      gc <- rbind(gc, data.frame(Gene = nm, Symbols = sy, stringsAsFactors = FALSE))
    }
    gene_catalog(gc)
    last_gene(nm)
    removeModal()
    showGenotypeModal(pref_gene = nm)
    showNotification(paste0("✅ Gene added: ", nm), type="message")
  })
  
  observeEvent(input$add_genotype, {
    req(current_cage())
    if (nrow(gene_catalog()) == 0) {
      # No genes configured -> open quick add
      showQuickAddGeneModal()
    } else {
      showGenotypeModal()
    }
  })
  
  observeEvent(input$move_cage, {
    if (!is.null(current_cage())) {
      showModal(modalDialog(
        title = paste("Move Mice from Cage:", current_cage()),
        div(
          textInput("new_cage_id", "New Cage ID:"),
          selectInput("move_reason", "Reason for Move:", 
                     choices = c("Breeding setup", "Experiment", "Health check", "Cage cleaning", "Other")),
          dateInput("move_date", "Move Date:", value = Sys.Date()),
          textInput("move_notes", "Notes:")
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("save_move", "Move Mice", class = "action-btn btn-edit")
        )
      ))
    }
  })
  
  # Save handlers for quick actions
  
  # Gene symbol UI output
  output$gene_symbol_ui <- renderUI({
    gene <- input$gene_choice %||% ""
    syms <- get_gene_syms(gene)
    if (!length(syms)) {
      return(textInput("genotype_result_manual",
                       "Genotype symbol (no symbols configured for this gene):", ""))
    } else {
      selectInput("genotype_result", "Genotype symbol:", choices = syms, selected = syms[1])
    }
  })
  
  # Remember last used gene
  observeEvent(input$gene_choice, {
    if (nzchar(input$gene_choice)) last_gene(input$gene_choice)
  })
  
  # Go to Genes tab quickly
  observeEvent(input$goto_genes_tab, {
    updateTabItems(session, "sidebar", "genes")
  })
  
  # New genotype save handler
  observeEvent(input$save_genotype_new, {
    req(current_cage())
    gene <- input$gene_choice %||% ""
    sym  <- input$genotype_result %||% input$genotype_result_manual %||% ""

    if (!nzchar(gene)) {           # ← NEW: no gene selected -> quick add dialog
      showQuickAddGeneModal(); return()
    }
    if (!nzchar(sym)) {
      showNotification("Pick a genotype symbol or type one.", type="warning"); return()
    }

    all_m <- mice_data_reactive()

    # Rows in current cage (absolute indices into all_m)
    in_cage_idx <- which(all_m$Cage == current_cage())

    scope <- input$apply_scope %||% "sel"

    if (scope == "cage") {
      cg <- cage_data_reactive()
      i  <- which(cg$Cage_Name == current_cage())[1]
      if (!is.na(i)) {
        cg$Genotype[i] <- set_gene_in_string(cg$Genotype[i], gene, sym)
        cg$Genotyping_Done[i] <- TRUE
        cg$Genotype_Confirmed[i] <- TRUE
        cage_data_reactive(cg)
      }
    } else {
      # NEW: map by stable IDs
      ids <- input$selected_mouse_ids %||% character(0)
      target_rows <- if (identical(scope, "all")) in_cage_idx else which(all_m$Mouse_ID %in% ids)
      if (!length(target_rows)) { 
        showNotification("No mice selected in this cage.", type="warning"); return() 
      }
      if (!"Genotype" %in% names(all_m)) all_m$Genotype <- ""
      for (r in target_rows) all_m$Genotype[r] <- set_gene_in_string(all_m$Genotype[r], gene, sym)
      mice_data_reactive(all_m)
    }
    
    # Log the action
    targets <- if (identical(input$apply_scope, "all")) sum(all_m$Cage == current_cage())
               else length(input$selected_mouse_ids %||% character(0))
    
    log_event(
      "Genotype Update",
      cage = current_cage(),
      count = targets,
      details = paste(input$gene_choice, "(", input$genotype_result %||% input$genotype_result_manual, ")",
                      "| scope:", input$apply_scope)
    )

    removeModal()
    showNotification(paste0("✅ Genotype set: ", gene, " (", sym, ")"), type="message")
  })
  
  # Legacy genotype save handler (keep for backward compatibility)
  observeEvent(input$save_genotype, {
    req(current_cage())
    cg <- cage_data_reactive()
    i <- which(cg$Cage_Name == current_cage())[1]
    if (!is.na(i)) {
      cg$Genotype[i] <- input$genotype_result
      cg$Genotyping_Done[i] <- TRUE
      cg$Genotype_Confirmed[i] <- isTRUE(input$genotype_confirmed)
      # Ensure Notes column exists
      if (!"Notes" %in% names(cg)) {
        cg$Notes <- rep("", nrow(cg))
      }
      add_note <- paste0(
        "Genotyped ", format(input$genotype_date, "%Y-%m-%d"),
        if (nzchar(input$genotype_lab)) paste0(" @", input$genotype_lab) else "",
        if (nzchar(input$genotype_notes)) paste0(" (", input$genotype_notes, ")") else ""
      )
      cg$Notes[i] <- gsub("^\\s+|\\s+$", "", paste(cg$Notes[i], add_note, sep = " | "))
      cage_data_reactive(cg)
    }
    removeModal()
    showNotification("✅ Genotype information saved!", type = "message")
  })
  
  observeEvent(input$save_move, {
    req(current_cage(), nzchar(input$new_cage_id))
    from <- current_cage()
    to <- input$new_cage_id
    all <- mice_data_reactive()
    if (sum(all$Cage == from) == 0) {
      showNotification("No mice in this cage to move.", type = "warning")
      return()
    }
    add_or_update_cage(to, status = "Holding", parent = from, setup = format(input$move_date, "%m/%d"))
    all$Cage[all$Cage == from] <- to
    mice_data_reactive(all)
    recalc_cage_counts()
    set_cage_sex_from_mice(c(from, to))
    
    # Log the action
    log_event(
      "Move Mice",
      cage = from,
      count = sum(all$Cage == to),   # after move it's the dest count
      details = paste("to", to, "| date:", format(input$move_date, "%Y-%m-%d"))
    )
    
    removeModal()
    showNotification(paste0("✅ Moved mice from ", from, " to ", to, "."), type = "message")
  })
  
  # Delete cage functionality
  observeEvent(input$delete_cage, {
    showModal(modalDialog(
      title = "⚠️ Confirm Deletion",
      div(
        p("Are you sure you want to delete this cage and all its mice?"),
        p("Cage:", current_cage()),
        p("This action cannot be undone.")
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_delete_cage", "Yes, Delete", class = "action-btn btn-danger")
      )
    ))
  })
  
  observeEvent(input$confirm_delete_cage, {
    current_cages <- cage_data_reactive()
    current_mice  <- mice_data_reactive()
    cage_to_delete <- current_cage()

    removed_count <- sum(current_mice$Cage == cage_to_delete)
    log_event("Delete Cage", cage = cage_to_delete, count = removed_count)

    current_mice  <- current_mice[current_mice$Cage != cage_to_delete, , drop = FALSE]
    current_cages <- current_cages[current_cages$Cage_Name != cage_to_delete, , drop = FALSE]

    mice_data_reactive(current_mice)
    cage_data_reactive(current_cages)
    recalc_cage_counts()
    current_cage(NULL)
    removeModal()
    showNotification(paste0("🗑️ Deleted cage ", cage_to_delete, " and ", removed_count, " mouse/mice."), type = "message")
  })
  
  # Delete mouse functionality
  observeEvent(input$delete_mouse, {
    showModal(modalDialog(
      title = "⚠️ Confirm Mouse Deletion",
      div(
        p("Are you sure you want to delete this mouse?"),
        p("Mouse ID:", current_mouse()),
        p("This action cannot be undone.")
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_delete_mouse", "Yes, Delete", class = "action-btn btn-danger")
      )
    ))
  })
  
  observeEvent(input$confirm_delete_mouse, {
    # Get current data
    current_mice <- mice_data_reactive()
    mouse_to_delete <- current_mouse()
    
    # Log the action before clearing state
    log_event(
      "Delete Mouse",
      mouse = mouse_to_delete,
      details = "Mouse record removed"
    )
    
    # Remove the mouse
    updated_mice <- current_mice[current_mice$Mouse_ID != mouse_to_delete, ]
    
    # Update reactive data
    mice_data_reactive(updated_mice)
    
    # Reset mouse selection
    current_mouse(NULL)
    
    removeModal()
    showNotification("✅ Mouse deleted successfully!", type = "message")
  })
  
  # Commit import logic
  observeEvent(input$commit_import, {
    df <- uploaded_df(); req(df)
    iss <- upload_issues()
    if (length(iss$errors)) {
      showNotification("❌ Fix errors before import.", type = "error")
      return()
    }

    # snapshot for undo
    prev_cages(cage_data_reactive())
    prev_mice(mice_data_reactive())

    if (identical(input$import_kind, "mice")) {
      # Determine an ID column present in uploaded df
      id_col <- c("Mouse_ID","Animal_ID","Physical Tag","Earmark")
      id_col <- id_col[id_col %in% names(df)][1]
      if (is.null(id_col)) {
        showNotification("❌ No mouse ID column found.", type = "error"); return()
      }

      cur <- mice_data_reactive()
      if (input$import_upsert) {
        # upsert by id_col -> Mouse_ID (if Mouse_ID not present, we add it as that id_col name)
        if (!"Mouse_ID" %in% names(df)) df$Mouse_ID <- df[[id_col]]

        # align columns: union of both sets
        all_cols <- union(names(cur), names(df))
        for (nm in setdiff(all_cols, names(cur))) cur[[nm]] <- ""
        for (nm in setdiff(all_cols, names(df)))  df[[nm]]  <- ""

        # index by Mouse_ID
        cur_ids <- cur$Mouse_ID
        up_ids  <- df$Mouse_ID

        # update existing rows (non-empty fields overwrite)
        in_both <- intersect(cur_ids, up_ids)
        for (id in in_both) {
          i <- which(cur_ids == id)
          j <- which(up_ids == id)
          for (nm in names(cur)) {
            val <- df[j, nm, drop=TRUE]
            if (!is.na(val) && nzchar(as.character(val))) cur[i, nm] <- val
          }
        }
        # append new rows
        to_add <- setdiff(up_ids, cur_ids)
        cur <- rbind(cur, df[df$Mouse_ID %in% to_add, names(cur), drop=FALSE])
        mice_data_reactive(cur)
      } else {
        # append only (no updates)
        if (!"Mouse_ID" %in% names(df)) df$Mouse_ID <- df[[id_col]]
        new_only <- df[!(df$Mouse_ID %in% mice_data_reactive()$Mouse_ID), , drop = FALSE]
        mice_data_reactive(bind_rows(mice_data_reactive(), new_only))
      }

      # Optionally recalc cage counts from mice
      if (isTRUE(input$import_recount)) {
        new_m <- mice_data_reactive()
        counts <- as.data.frame(table(new_m$Cage), stringsAsFactors = FALSE)
        names(counts) <- c("Cage_Name","Num_Mice")
        cg <- cage_data_reactive()
        cg <- merge(cg[, setdiff(names(cg), "Num_Mice"), drop=FALSE], counts, by = "Cage_Name", all.x = TRUE)
        cg$Num_Mice[is.na(cg$Num_Mice)] <- 0
        cage_data_reactive(cg)
      }

      # Log the import
      log_event("Import Mice",
                details = paste("Rows:", nrow(df), "| Upsert:", isTRUE(input$import_upsert)))
      
      showNotification("✅ Mice imported.", type = "message")
      removeModal()
    } else {
      # cage-based import: upsert by Cage_Name
      df$Cage_Name <- as.character(df$Cage_Name)
      cur <- cage_data_reactive()

      # align columns
      all_cols <- union(names(cur), names(df))
      for (nm in setdiff(all_cols, names(cur))) cur[[nm]] <- ""
      for (nm in setdiff(all_cols, names(df)))  df[[nm]]  <- ""

      if (input$import_upsert) {
        cur_names <- cur$Cage_Name
        up_names  <- df$Cage_Name
        in_both <- intersect(cur_names, up_names)
        for (id in in_both) {
          i <- which(cur_names == id)
          j <- which(up_names == id)
          for (nm in names(cur)) {
            val <- df[j, nm, drop=TRUE]
            if (!is.na(val) && nzchar(as.character(val))) cur[i, nm] <- val
          }
        }
        to_add <- setdiff(up_names, cur_names)
        cur <- rbind(cur, df[df$Cage_Name %in% to_add, names(cur), drop=FALSE])
        cage_data_reactive(cur)
      } else {
        to_add <- df[!(df$Cage_Name %in% cur$Cage_Name), names(cur), drop=FALSE]
        cage_data_reactive(bind_rows(cur, to_add))
      }
      # Log the import
      log_event("Import Cages",
                details = paste("Rows:", nrow(df), "| Upsert:", isTRUE(input$import_upsert)))
      
      showNotification("✅ Cages imported.", type = "message")
      removeModal()
    }
  })
  
  # Undo last import (optional)
  observeEvent(input$undo_import, {
    if (!is.null(prev_cages())) cage_data_reactive(prev_cages())
    if (!is.null(prev_mice()))  mice_data_reactive(prev_mice())
    showNotification("↩️ Reverted to pre-import state.", type="message")
  })
  
  # Back to list controls for Cage View
  observeEvent(input$back_to_list_btn,  { current_cage(NULL) })
  observeEvent(input$back_to_cage_list, { current_cage(NULL) })  # from breadcrumb click
  
  # Bulk delete handlers moved below with confirmation modals
  
  # ----- Bulk delete: MICE (Mouse Records page) -----
  observeEvent(input$delete_selected_mice_main, {
    ids <- input$mouse_selected_keys_main %||% character(0)
    if (!length(ids)) {
      showNotification("Select one or more mice using the left checkboxes.", type="warning"); return()
    }
    showModal(modalDialog(
      title = "Delete selected mice?",
      HTML(paste("IDs:", paste(ids, collapse = ", "))),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_delete_mice_main", "Delete", class="action-btn btn-danger")
      )
    ))
  })

  observeEvent(input$confirm_delete_mice_main, {
    ids <- input$mouse_selected_keys_main %||% character(0)
    if (!length(ids)) return()
    mice <- mice_data_reactive()
    for (m in ids) log_event("Delete Mouse", mouse = m, details = "Bulk delete")
    mice <- mice[!(mice$Mouse_ID %in% ids), , drop = FALSE]
    mice_data_reactive(mice)
    recalc_cage_counts()
    removeModal()
    showNotification(paste0("✅ Deleted ", length(ids), " mouse/mice."), type="message")
  })
  
  # --- Bulk delete: selected mice inside current cage --------------------------
  observeEvent(input$delete_selected_mice_in_cage, {
    ids <- input$selected_mouse_ids %||% character(0)
    if (!length(ids)) { showNotification("Select one or more mice first.", type="warning"); return() }

    all <- mice_data_reactive()
    mice_data_reactive(all[!(all$Mouse_ID %in% ids), , drop=FALSE])
    recalc_cage_counts()
    set_cage_sex_from_mice(current_cage())

    log_event("Delete Mice (Cage)", cage = current_cage(),
              count = length(ids), details = paste(ids, collapse = ", "))
    showNotification(paste0("🗑️ Deleted ", length(ids), " mouse record(s)."), type="message")
  })
  
  # ----- Bulk delete: CAGES (Home) -----
  observeEvent(input$delete_selected_cages_dashboard, {
    keys <- input$cage_selected_keys_dashboard %||% character(0)
    if (!length(keys)) {
      showNotification("Select one or more cages using the left checkboxes.", type="warning"); return()
    }
    showModal(modalDialog(
      title = "Delete selected cages?",
      HTML(paste("Cages:", paste(keys, collapse = ", "))),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_delete_cages_dashboard", "Delete", class="action-btn btn-danger")
      )
    ))
  })

  observeEvent(input$confirm_delete_cages_dashboard, {
    keys <- input$cage_selected_keys_dashboard %||% character(0)
    if (!length(keys)) return()
    cages <- cage_data_reactive()
    mice  <- mice_data_reactive()

    for (k in keys) log_event("Delete Cage", cage = k, count = sum(mice$Cage == k), details = "Bulk delete")

    mice  <- mice[!(mice$Cage %in% keys), , drop = FALSE]
    cages <- cages[!(cages$Cage_Name %in% keys), , drop = FALSE]

    cage_data_reactive(cages)
    mice_data_reactive(mice)
    recalc_cage_counts()
    current_cage(NULL)
    removeModal()
    showNotification(paste0("✅ Deleted ", length(keys), " cage(s)."), type="message")
  })

  # ----- Bulk delete: CAGES (Cage View list) -----
  observeEvent(input$delete_selected_cages_cv, {
    keys <- input$cage_selected_keys_cv %||% character(0)
    if (!length(keys)) {
      showNotification("Select one or more cages using the left checkboxes.", type="warning"); return()
    }
    showModal(modalDialog(
      title = "Delete selected cages?",
      HTML(paste("Cages:", paste(keys, collapse = ", "))),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_delete_cages_cv", "Delete", class="action-btn btn-danger")
      )
    ))
  })

  observeEvent(input$confirm_delete_cages_cv, {
    keys <- input$cage_selected_keys_cv %||% character(0)
    if (!length(keys)) return()
    cages <- cage_data_reactive()
    mice  <- mice_data_reactive()

    for (k in keys) log_event("Delete Cage", cage = k, count = sum(mice$Cage == k), details = "Bulk delete")

    mice  <- mice[!(mice$Cage %in% keys), , drop = FALSE]
    cages <- cages[!(cages$Cage_Name %in% keys), , drop = FALSE]

    cage_data_reactive(cages)
    mice_data_reactive(mice)
    recalc_cage_counts()
    current_cage(NULL)
    removeModal()
    showNotification(paste0("✅ Deleted ", length(keys), " cage(s)."), type="message")
  })
  
  # Handle various button clicks
  observeEvent(input$add_cage, {
    showModal(modalDialog(
      title = "Add New Cage",
      div(
        textInput("new_cage_id", "Cage ID:"),
        selectInput("new_cage_genotype", "Genotype:", 
                   choices = c("Cx3cr1", "Tmem119", "Other")),
        selectInput("new_cage_sex", "Sex:", choices = c("M", "F", "Mixed")),
        textInput("new_cage_notes", "Notes:")
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_new_cage", "Save", class = "action-btn btn-view")
      )
    ))
  })
  
  observeEvent(input$export_cages, {
    # Simulate export functionality
    showNotification("📤 Exporting cage data to Excel...", type = "message")
  })
  
  observeEvent(input$generate_report, {
    # Simulate report generation
    showNotification("📊 Generating colony report...", type = "message")
  })
  
  observeEvent(input$print_view, {
    # Simulate print functionality
    showNotification("🖨️ Opening print view...", type = "message")
  })
  
  # Add loading states
  observeEvent(input$save_new_cage, {
    req(input$new_cage_id)
    add_or_update_cage(
      cage_name = input$new_cage_id,
      genotype  = input$new_cage_genotype %||% "Unknown",
      sex       = input$new_cage_sex %||% "Mixed",
      status    = input$new_cage_status %||% "Active",
      parent    = "",
      setup     = format(Sys.Date(), "%m/%d")
    )
    recalc_cage_counts()
    removeModal()
    showNotification("✅ Cage added successfully!", type = "message")
  })
  
  # Open the Upload / Import modal from the Mouse Records page
  observeEvent(input$import_mice, {
    showModal(uploadModalUI(!is.null(mlims_cols())))
  })
  
  # Show Excel-specific options after a file is chosen (if Excel)
  observeEvent(input$upload_file, {
    req(input$upload_file$datapath)
    ext <- tolower(tools::file_ext(input$upload_file$name))
    if (ext %in% c("xlsx","xls")) {
      # Grab sheet names
      sheets <- readxl::excel_sheets(input$upload_file$datapath)
      output$excel_options_ui <- renderUI({
        tagList(
          selectInput("excel_sheet", "Sheet", choices = sheets, selected = sheets[1]),
          numericInput("excel_skip", "Rows to skip before header", value = 0, min = 0, max = 10000),
          checkboxInput("excel_has_header", "First non-skipped row contains column names", TRUE)
        )
      })
    } else {
      output$excel_options_ui <- renderUI(NULL)
    }
  })
  
  # Clear paste box
  observeEvent(input$clear_paste, {
    updateTextAreaInput(session, "paste_table", value = "")
  })
  
  # Handle pasted table -> unify with the same preview/validation
  observeEvent(input$paste_table, {
    txt <- input$paste_table
    if (!nzchar(txt)) {
      uploaded_df(NULL)
      upload_issues(list(errors = character(), warnings = character()))
      return()
    }

    df <- read_pasted_df(txt)
    if (is.null(df) || !nrow(df)) {
      uploaded_df(NULL)
      upload_issues(list(errors = "Could not parse pasted table. Ensure the first row is headers.", warnings = character()))
      return()
    }

    # Shape & validate like the file path flow
    parse_and_stage_upload(df)
  })
  
  # Show Excel-specific options after a file is chosen (if Excel)
  observeEvent(input$upload_file, {
    req(input$upload_file$datapath)
    ext <- tolower(tools::file_ext(input$upload_file$name))
    if (ext %in% c("xlsx","xls")) {
      # Grab sheet names
      sheets <- readxl::excel_sheets(input$upload_file$datapath)
      output$excel_options_ui <- renderUI({
        tagList(
          selectInput("excel_sheet", "Sheet", choices = sheets, selected = sheets[1]),
          numericInput("excel_skip", "Rows to skip before header", value = 0, min = 0, max = 10000),
          checkboxInput("excel_has_header", "First non-skipped row contains column names", TRUE)
        )
      })
    } else {
      output$excel_options_ui <- renderUI(NULL)
    }
  })
  
  # Clear paste box
  observeEvent(input$clear_paste, {
    updateTextAreaInput(session, "paste_table", value = "")
  })
  
  # Unified "parse upload" path with progress
  parse_and_stage_upload <- function(df) {
    if (identical(input$import_kind, "mice")) {
      # Shape to mLIMS headers if present
      if (!is.null(mlims_cols())) {
        cols <- mlims_cols()
        out <- as.data.frame(
          setNames(replicate(length(cols), rep("", nrow(df)), simplify = FALSE), cols),
          stringsAsFactors = FALSE
        )
        common <- intersect(names(df), cols)
        for (nm in common) out[[nm]] <- as.character(df[[nm]])
        df <- out
      }
      v <- validate_mice_df(df)
      upload_issues(v)
      uploaded_df(df)

    } else { # cage-based
      keep <- intersect(names(df), cage_template_cols)
      df2 <- as.data.frame(df[keep], stringsAsFactors = FALSE, check.names = FALSE)
      missing <- setdiff(cage_template_cols, names(df2))
      for (m in missing) df2[[m]] <- ""
      df2 <- df2[, cage_template_cols, drop = FALSE]
      v <- validate_cage_df(df2)
      upload_issues(v)
      uploaded_df(df2)
    }
  }
  
  # Read the chosen file and stage it for preview/import
  observeEvent(input$upload_file, {
    req(input$upload_file$datapath)
    withProgress(message = "Reading file…", value = 0, {
      incProgress(0.2)
      ext <- tolower(tools::file_ext(input$upload_file$name))

      if (ext %in% c("xlsx","xls")) {
        # Excel: use options if present, otherwise defaults
        sheet <- input$excel_sheet %||% 1
        skip  <- input$excel_skip %||% 0
        hdr   <- if (is.null(input$excel_has_header)) TRUE else isTRUE(input$excel_has_header)
        df <- read_upload_df(input$upload_file$datapath, sheet = sheet, skip = skip, col_names = hdr)
      } else {
        df <- read_upload_df(input$upload_file$datapath)
      }
      incProgress(0.6)
      if (is.null(df) || !nrow(df)) {
        uploaded_df(NULL)
        upload_issues(list(errors = "No rows found in the uploaded file.", warnings = character()))
        return()
      }
      parse_and_stage_upload(df)
      incProgress(0.2)
    })
  })
  
  # If user changes Excel options after upload, re-parse automatically
  observeEvent(list(input$excel_sheet, input$excel_skip, input$excel_has_header, input$import_kind), {
    req(input$upload_file$datapath)
    ext <- tolower(tools::file_ext(input$upload_file$name))
    if (!(ext %in% c("xlsx","xls"))) return()

    withProgress(message = "Applying Excel options…", value = 0, {
      sheet <- input$excel_sheet %||% 1
      skip  <- input$excel_skip %||% 0
      hdr   <- if (is.null(input$excel_has_header)) TRUE else isTRUE(input$excel_has_header)
      df <- read_upload_df(input$upload_file$datapath, sheet = sheet, skip = skip, col_names = hdr)
      incProgress(0.6)
      if (is.null(df) || !nrow(df)) {
        uploaded_df(NULL)
        upload_issues(list(errors = "No rows found after applying Excel options.", warnings = character()))
        return()
      }
      parse_and_stage_upload(df)
      incProgress(0.4)
    })
  })
  
  # Handle pasted table -> unify with the same preview/validation
  observeEvent(input$paste_table, {
    txt <- input$paste_table
    if (!nzchar(txt)) {
      uploaded_df(NULL)
      upload_issues(list(errors = character(), warnings = character()))
      return()
    }

    df <- read_pasted_df(txt)
    if (is.null(df) || !nrow(df)) {
      uploaded_df(NULL)
      upload_issues(list(errors = "Could not parse pasted table. Ensure the first row is headers.", warnings = character()))
      return()
    }

    # Shape & validate like the file path flow
    parse_and_stage_upload(df)
  })
  
  # Fast preview for big data
  output$import_preview <- DT::renderDT({
    df <- uploaded_df()
    if (is.null(df)) {
      return(DT::datatable(data.frame(Upload = "No file yet"), options = list(dom='t')))
    }
    n <- input$preview_rows %||% 200
    head_df <- utils::head(df, n)
    DT::datatable(
      head_df,
      extensions = "Buttons",
      options = list(pageLength = 10, dom = 'Bfrtip',
                     buttons = c('copy','csv','excel','pdf','print'),
                     server = TRUE),
      rownames = FALSE
    )
  })
  
  # --- Daily Log server components ---
  
  # Action filter choices always reflect current log
  output$log_action_filter_ui <- renderUI({
    df <- audit_log()
    acts <- sort(unique(df$Action))
    selectizeInput("log_action_filter", "Action filter (optional)", choices = c("All", acts), selected = "All")
  })

  daily_log_for_date <- reactive({
    df <- audit_log()
    if (!nrow(df)) return(df)
    d  <- as.Date(input$log_date %||% Sys.Date())
    keep <- as.Date(df$Time) == d
    df <- df[keep, , drop = FALSE]
    if (!is.null(input$log_action_filter) && input$log_action_filter != "All") {
      df <- df[df$Action == input$log_action_filter, , drop = FALSE]
    }
    df[order(df$Time), , drop = FALSE]
  })

  output$daily_log_table <- DT::renderDT({
    df <- daily_log_for_date()
    if (!nrow(df)) {
      return(DT::datatable(data.frame(Message = "No activity logged for this date."),
                           options = list(dom='t'), rownames = FALSE))
    }
    pretty <- transform(df,
                        Time = format(Time, "%Y-%m-%d %H:%M"),
                        Count = ifelse(is.na(Count), "", as.character(Count)))
    DT::datatable(pretty,
      extensions = c("Buttons"),
      options = list(dom = 'Bfrtip', buttons = c('copy','csv','excel','pdf','print'),
                     pageLength = 10, order = list(list(0, 'asc'))),
      rownames = FALSE
    )
  })

  output$daily_log_summary_ui <- renderUI({
    df <- daily_log_for_date()
    if (!nrow(df)) return(NULL)
    sum_by <- aggregate(rep(1, nrow(df)), by = list(Action = df$Action), FUN = sum)
    htmltools::HTML(paste0(
      "<h4>Summary (", format(input$log_date, "%Y-%m-%d"), ")</h4>",
      "<ul>",
      paste(apply(sum_by, 1, function(r) sprintf("<li>%s: %s</li>", htmltools::htmlEscape(r[1]), r[2])), collapse = ""),
      "</ul>"
    ))
  })

  # Printable HTML (self-contained)
  create_daily_log_html <- function(df, date_label) {
    esc <- function(z) htmltools::htmlEscape(as.character(z %||% ""))
    rows <- if (nrow(df)) paste(apply(df, 1, function(r) {
      sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td style='text-align:right;'>%s</td><td>%s</td></tr>",
              esc(format(r[["Time"]], "%Y-%m-%d %H:%M")),
              esc(r[["Action"]]), esc(r[["Cage"]]), esc(r[["Mouse_ID"]]),
              esc(ifelse(is.na(r[["Count"]]), "", r[["Count"]])),
              esc(r[["Details"]]))
    }), collapse = "") else ""

    # summary
    sum_by <- if (nrow(df)) aggregate(rep(1, nrow(df)), by = list(Action = df$Action), FUN = sum) else NULL
    sum_html <- if (!is.null(sum_by) && nrow(sum_by)) {
      paste0("<ul>", paste(apply(sum_by,1, function(r) sprintf("<li>%s: %s</li>", esc(r[1]), r[2])), collapse=""), "</ul>")
    } else "<p>No activity recorded for this day.</p>"

    paste0(
      sprintf("<h1>mLIMS — Daily Log</h1><div class='meta'>Date: %s</div>", esc(date_label)),
      "<table><thead><tr><th>Time</th><th>Action</th><th>Cage</th><th>Mouse</th><th>Count</th><th>Details</th></tr></thead><tbody>",
      paste(rows, collapse = "\n"),
      "</tbody></table>",
      "</body></html>"
    )
  }

  output$dl_daily_log_html <- downloadHandler(
    filename = function() sprintf("mLIMS_daily_log_%s.html", format(input$log_date %||% Sys.Date(), "%Y-%m-%d")),
    content = function(file) {
      df <- daily_log_for_date()
      html <- create_daily_log_html(df, format(input$log_date %||% Sys.Date(), "%Y-%m-%d"))
      writeLines(html, con = file, useBytes = TRUE)
    }
  )

  # Full raw audit CSV (all days)
  output$dl_audit_csv <- downloadHandler(
    filename = function() sprintf("mLIMS_audit_%s.csv", format(Sys.Date(), "%Y-%m-%d")),
    content = function(file) {
      df <- audit_log()
      df$Time <- format(df$Time, "%Y-%m-%d %H:%M:%S")
      utils::write.csv(df, file, row.names = FALSE, na = "")
    }
  )
  
  # --- Database Persistence (Debounced) ---
  
  # Debounced persist: saves after the user stops changing things for ~1s
  cages_for_save <- reactive({ cage_data_reactive() })
  mice_for_save  <- reactive({ mice_data_reactive() })
  genes_for_save <- reactive({ gene_catalog() })
  audit_for_save <- reactive({ audit_log() })

  cages_deb <- debounce(cages_for_save, 1000)
  mice_deb  <- debounce(mice_for_save,  1000)
  genes_deb <- debounce(genes_for_save, 1000)

  observeEvent(list(cages_deb(), mice_deb(), genes_deb()), {
    # Do not include the whole audit log each time—only append NEW rows
    # Here we append nothing (audit handled below), just replace 3 main tables.
    persist_all(pool,
                cages = cages_deb(),
                mice  = mice_deb(),
                genes = genes_deb(),
                audit_append = NULL)
  }, ignoreInit = FALSE)

  # Append audit events as they arrive (not replace the whole table)
  last_audit_rows <- reactiveVal(0L)
  observeEvent(audit_for_save(), {
    df <- audit_for_save()
    n  <- nrow(df)
    if (n > last_audit_rows()) {
      new <- df[seq.int(last_audit_rows() + 1L, n), , drop = FALSE]
      persist_all(pool, cages = NULL, mice = NULL, genes = NULL, audit_append = new)
      last_audit_rows(n)
    }
  }, ignoreInit = FALSE)
  
  # --- Backup & Restore Handlers ---
  
  observeEvent(input$backup_now, {
    if (is.null(s3_board)) { showNotification("S3 backup not configured.", type = "error"); return() }
    payload <- list(
      cages = cage_data_reactive(),
      mice  = mice_data_reactive(),
      genes = gene_catalog(),
      audit = audit_log(),
      ts    = Sys.time()
    )
    pins::pin_write(s3_board, payload, name = "mlims_backup", type = "rds")
    showNotification("☁️ Backup saved to S3 (versioned).", type = "message")
  })

  observeEvent(input$restore_latest, {
    if (is.null(s3_board)) { showNotification("S3 backup not configured.", type = "error"); return() }
    payload <- pins::pin_read(s3_board, "mlims_backup")
    if (!is.null(payload$cages)) cage_data_reactive(payload$cages)
    if (!is.null(payload$mice))  mice_data_reactive(payload$mice)
    if (!is.null(payload$genes)) gene_catalog(payload$genes)
    if (!is.null(payload$audit)) audit_log(payload$audit)
    recalc_cage_counts()
    showNotification("🔁 Restored latest backup.", type = "message")
  })
  
  # Database status UI
  output$db_status_ui <- renderUI({
    if (is.null(pool)) return(div(style = "color:#b71c1c;", "Database not configured."))
    ok <- tryCatch({ DBI::dbGetQuery(pool, "SELECT 1 AS ok")$ok[1] == 1 }, error = function(e) FALSE)
    if (ok) {
      tagList(
        div(style = "color:#2e7d32;", HTML("✅ Connected to Postgres.")),
        tags$small(style = "color:#555;", "Edits are saved to the DB and included in backups.")
      )
    } else {
      div(style = "color:#ef6c00;", HTML("⚠️ Database reachable but ping failed."))
    }
  })
  
  # Heartbeat to keep session alive
  output$heartbeat <- renderText({
    invalidateLater(30000, session)  # ping every 30s
    ""                               # no visible output
  })
  
  # Initialize database connection and load data (moved here)
  db_state <- load_state_from_db()
  if (!is_db_empty(db_state$cages))  cage_data_reactive(db_state$cages)
  if (!is_db_empty(db_state$mice))   mice_data_reactive(db_state$mice)
  if (!is_db_empty(db_state$genes))  gene_catalog(db_state$genes)
  
  # Save state whenever tab or cage changes
  observeEvent(input$sidebar, {
    session$sendCustomMessage("rememberState", list(
      tab  = input$sidebar,
      cage = current_cage()
    ))
  }, ignoreInit = FALSE)

  observeEvent(current_cage(), {
    session$sendCustomMessage("rememberState", list(
      tab  = input$sidebar,
      cage = current_cage()
    ))
  }, ignoreInit = FALSE)

  # Restore handlers (these inputs are set by the JS above)
  observeEvent(input$restore_tab, {
    if (!is.null(input$restore_tab)) {
      updateTabItems(session, "sidebar", input$restore_tab)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$restore_cage, {
    if (!is.null(input$restore_cage) && nzchar(input$restore_cage)) {
      current_cage(input$restore_cage)
      updateTabItems(session, "sidebar", "cage_view")
    }
  }, ignoreInit = TRUE)
  
  # Navigate breadcrumb back to mice tab if clicked in modal
  observeEvent(input$breadcrumb_click, {
    if (identical(input$breadcrumb_click, "mice")) updateTabItems(session, "sidebar", "mice")
  }, ignoreInit = TRUE)
}

# Run the application
shinyApp(ui = ui, server = server)
