################################################################################
# qPCR Cleaner and Analyzer
# Enhanced version with improved UI, error handling, and user experience
#
# Features:
# - Upload CSV or Excel files with or without headers
# - Automatic data validation and cleaning
# - Flexible replicate handling (3 or 4 replicates)
# - Outlier detection and removal
# - ΔCt calculation with optimal pairing
# - Expression analysis and relative expression
# - Statistical summaries with mean and SD
# - Interactive data visualization
# - Comprehensive error handling
# - Clear user instructions
################################################################################

library(shiny)
library(gtools)
library(readxl)
library(DT)
library(dplyr)
library(tidyr)
library(shinyjs)
library(shinyWidgets)

# Custom CSS for better styling
custom_css <- "
.nav-tabs > li > a {
  background-color: #f8f9fa;
  border: 1px solid #dee2e6;
  margin-right: 2px;
}

.nav-tabs > li.active > a {
  background-color: #007bff;
  color: white;
  border: 1px solid #007bff;
}

.well {
  background-color: #f8f9fa;
  border: 1px solid #dee2e6;
  border-radius: 5px;
  padding: 15px;
  margin-bottom: 15px;
}

.btn-primary {
  background-color: #007bff;
  border-color: #007bff;
}

.btn-primary:hover {
  background-color: #0056b3;
  border-color: #0056b3;
}

.alert {
  border-radius: 5px;
  margin-bottom: 15px;
}

.table-responsive {
  border-radius: 5px;
  overflow: hidden;
}

.progress {
  height: 20px;
  border-radius: 10px;
}

.progress-bar {
  background-color: #007bff;
}
"

ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$style(HTML(custom_css))),
  
  tags$script(HTML("
    document.addEventListener('DOMContentLoaded', function(){
      Shiny.addCustomMessageHandler('clipboardRead', function(msg){
        if (navigator.clipboard && navigator.clipboard.readText) {
          navigator.clipboard.readText().then(function(txt){
            Shiny.setInputValue('clip_payload', txt, {priority: 'event'});
          }).catch(function(e){
            Shiny.setInputValue('clip_payload', 'ERROR: ' + e, {priority: 'event'});
          });
        } else {
          Shiny.setInputValue('clip_payload', 'ERROR: Clipboard API not available', {priority: 'event'});
        }
      });
    });
  ")),
  
  titlePanel(
    div(
      h1("qPCR Cleaner and Analyzer", style = "color: #2c3e50; margin-bottom: 5px;"),
      h4("Process and analyze qPCR data with advanced statistical methods", 
         style = "color: #7f8c8d; font-weight: normal; margin-top: 0;")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      wellPanel(
        h4("📁 Data Upload", style = "color: #2c3e50; margin-top: 0;"),
        
        fileInput("file", 
                 label = "Upload your qPCR data file:",
                 accept = c(".csv", ".xlsx", ".xls"),
                 placeholder = "Choose CSV or Excel file..."),
        
        checkboxInput("has_header", "File has column headers", value = FALSE),
        
        tags$hr(),
        h4("📋 Or paste data"),
        helpText("Paste a numeric grid (Excel/CSV). Rows = samples; columns = gene blocks of replicates (ref gene first)."),
        div(
          class = "btn-group", style="width:100%; margin-bottom:8px;",
          actionButton("paste_clip", "Paste from clipboard", class = "btn btn-default"),
          checkboxInput("paste_transpose", "Transpose pasted grid", FALSE)
        ),
        textAreaInput(
          "paste_box", NULL, placeholder = "Paste here… (Ctrl/Cmd+V)",
          rows = 8, resize = "vertical"
        ),
        checkboxInput("paste_has_header", "First row is headers", FALSE),
        checkboxInput("paste_has_left_labels", "First column has sample labels", FALSE),
        actionButton("load_paste", "Use pasted data", class = "btn-primary", style = "width:100%;"),
        
        hr(),
        
        h5("🔧 Analysis Settings", style = "color: #2c3e50;"),
        
        radioButtons("reps", 
                    "Number of replicates per gene:",
                    choices = c("3 replicates" = 3, "4 replicates" = 4),
                    selected = 3,
                    inline = TRUE),
        
        checkboxInput("remove_outliers", 
                     "Remove outliers (for 4 replicates)", 
                     value = TRUE),
        
        numericInput("outlier_threshold", 
                    "Outlier detection threshold (SD):",
                    value = 2.0, 
                    min = 1.0, 
                    max = 5.0, 
                    step = 0.1),
        
        hr(),
        
        actionButton("analyze", 
                    "🚀 Run Analysis", 
                    class = "btn-primary btn-lg",
                    style = "width: 100%; margin-bottom: 10px;"),
        
        actionButton("reset", 
                    "🔄 Reset", 
                    class = "btn-secondary",
                    style = "width: 100%;"),
        
        hr(),
        
        conditionalPanel(
          condition = "output.analysis_complete",
          div(style="display:flex; gap:8px; width:100%;",
            downloadButton("download_csv_stacked", "📥 Download Results (stacked CSV)", class="btn-success", style="flex:1;"),
            downloadButton("download_csv_tidy",    "📥 Download Tidy CSV",             class="btn-default", style="flex:1;")
          )
        )
      )
    ),
    
    mainPanel(
      width = 9,
      
      tabsetPanel(
        id = "main_tabs",
        type = "tabs",
        
        # Instructions Tab
        tabPanel(
          "📖 Instructions",
          wellPanel(
            h3("How to Use This App", style = "color: #2c3e50;"),
            
            h4("📋 Data Format Requirements:"),
            tags$ul(
              tags$li("Upload CSV or Excel files containing your qPCR Ct values"),
              tags$li("Data should be organized with genes in columns and samples in rows"),
              tags$li("First gene should be your reference/control gene"),
              tags$li("Missing values can be marked as 'Undetermined', '#VALUE!', or left blank"),
              tags$li("Choose whether your file has column headers")
            ),
            
            h4("⚙️ Analysis Options:"),
            tags$ul(
              tags$li("Select number of replicates per gene (3 or 4)"),
              tags$li("For 4 replicates, optionally remove outliers automatically"),
              tags$li("Adjust outlier detection sensitivity if needed")
            ),
            
            h4("📊 Output Information:"),
            tags$ul(
              tags$li("Cleaned Ct values with outliers removed"),
              tags$li("ΔCt values calculated using optimal replicate pairing"),
              tags$li("Expression values (2^ΔCt)"),
              tags$li("Relative expression compared to control"),
              tags$li("Statistical summaries (mean ± SD)")
            ),
            
            hr(),
            
            h4("💡 Tips:"),
            tags$ul(
              tags$li("Ensure your reference gene is in the first column"),
              tags$li("Check the data preview before running analysis"),
              tags$li("Review outlier detection results"),
              tags$li("Download results for further analysis in other software")
            )
          )
        ),
        
        # Data Preview Tab
        tabPanel(
          "👀 Data Preview",
          wellPanel(
            h3("Raw Data Preview", style = "color: #2c3e50;"),
            
            conditionalPanel(
              condition = "!output.file_uploaded",
              div(
                class = "alert alert-info",
                icon("info-circle"),
                "Please upload a file to see the data preview."
              )
            ),
            
            conditionalPanel(
              condition = "output.file_uploaded",
              div(
                style = "margin-bottom: 15px;",
                verbatimTextOutput("file_info")
              ),
              
              div(
                style = "margin-bottom: 15px;",
                h5("First 10 rows of data:"),
                DTOutput("raw_data_table")
              ),
              
              conditionalPanel(
                condition = "output.has_outliers",
                div(
                  class = "alert alert-warning",
                  icon("exclamation-triangle"),
                  "Outliers detected in your data. Consider enabling outlier removal."
                )
              )
            )
          )
        ),
        
        # Analysis Results Tab
        tabPanel(
          "📈 Analysis Results",
          conditionalPanel(
            condition = "!output.analysis_complete",
            wellPanel(
              div(
                class = "alert alert-info",
                icon("info-circle"),
                "Click 'Run Analysis' to process your data and view results."
              )
            )
          ),
          
          conditionalPanel(
            condition = "output.analysis_complete",
            
            # Summary Statistics
            wellPanel(
              h3("📊 Summary Statistics", style = "color: #2c3e50;"),
              
              fluidRow(
                column(6,
                  h5("Expression Analysis:"),
                  tableOutput("expr_summary_table")
                ),
                column(6,
                  h5("Relative Expression Analysis:"),
                  tableOutput("rel_summary_table")
                )
              )
            ),
            
            # Detailed Results
            wellPanel(
              h3("🔍 Detailed Results", style = "color: #2c3e50;"),
              
              tabsetPanel(
                type = "pills",
                
                tabPanel("Cleaned Ct Values",
                  DTOutput("ct_table")
                ),
                
                tabPanel("ΔCt Values",
                  DTOutput("dct_table")
                ),
                
                tabPanel("Expression Values",
                  DTOutput("expr_table")
                ),
                
                tabPanel("Relative Expression",
                  DTOutput("rel_table")
                ),
                
                tabPanel("Relative (ratio of means)",
                  DTOutput("rel_means_table")
                )
              )
            ),
            

          )
        ),
        
        # Quality Control Tab
        tabPanel(
          "🔬 Quality Control",
          conditionalPanel(
            condition = "!output.analysis_complete",
            wellPanel(
              div(
                class = "alert alert-info",
                icon("info-circle"),
                "Run analysis to view quality control metrics."
              )
            )
          ),
          
          conditionalPanel(
            condition = "output.analysis_complete",
            
            wellPanel(
              h3("Quality Control Metrics", style = "color: #2c3e50;"),
              
              fluidRow(
                column(6,
                  h5("Data Quality Summary:"),
                  tableOutput("qc_summary")
                ),
                column(6,
                  h5("Outlier Detection Results:"),
                  tableOutput("outlier_summary")
                )
              ),
              

            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive values to store data
  rv <- reactiveValues(
    raw_data = NULL,
    cleaned_data = NULL,
    dct_data = NULL,
    expr_data = NULL,
    rel_data = NULL,
    csv_data = NULL,
    analysis_complete = FALSE,
    file_uploaded = FALSE,
    has_outliers = FALSE,
    qc_metrics = NULL,
    outlier_info = NULL
  )
  
  # Always use an integer for replicate count
  reps_num <- reactive({ as.integer(input$reps) })
  
  # Robust paste parser
  parse_pasted_matrix <- function(txt, has_header = FALSE, has_left_labels = FALSE, transpose = FALSE) {
    stopifnot(is.character(txt), length(txt) == 1)
    if (nchar(trimws(txt)) == 0) stop("No text to parse")

    # Normalize line endings
    txt <- gsub("\r", "", txt)

    # Guess separator: prefer tab, then comma, then semicolon, else whitespace
    sep <- if (grepl("\t", txt)) "\t" else if (grepl(",", txt)) "," else if (grepl(";", txt)) ";" else ""

    # Read table
    con <- textConnection(txt); on.exit(close(con), add = TRUE)
    df <- tryCatch({
      read.table(con, sep = sep, header = has_header, quote = "\"'", comment.char = "",
                 stringsAsFactors = FALSE, na.strings = c("Undetermined", "#VALUE!", "NA", ""))
    }, error = function(e) stop("Could not parse pasted text: ", e$message))

    # Drop all-NA rows/cols
    df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]
    df <- df[, colSums(!is.na(df)) > 0, drop = FALSE]

    # Optional: strip left label column and/or top header row if user indicated
    if (isTRUE(has_left_labels) && ncol(df) > 0) df <- df[, -1, drop = FALSE]

    # Heuristics: if top row is mostly non-numeric, drop it; if left col mostly non-numeric, drop it
    is_num_row <- function(v) mean(!is.na(suppressWarnings(as.numeric(v)))) >= 0.5
    is_num_col <- function(v) mean(!is.na(suppressWarnings(as.numeric(v)))) >= 0.5
    if (nrow(df) > 0 && !is_num_row(df[1, ])) df <- df[-1, , drop = FALSE]
    while (ncol(df) > 0 && !is_num_col(df[, 1])) df <- df[, -1, drop = FALSE]

    # Coerce to numeric
    df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))

    # Drop any all-NA rows/cols after coercion
    df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]
    df <- df[, colSums(!is.na(df)) > 0, drop = FALSE]

    if (isTRUE(transpose)) df <- as.data.frame(t(as.matrix(df)))

    if (nrow(df) == 0 || ncol(df) == 0) stop("No numeric cells detected after cleaning")
    df
  }
  
  # File upload observer
  observeEvent(input$file, {
    req(input$file)
    
    tryCatch({
      file_path <- input$file$datapath
      file_name <- input$file$name
      
      # Read file based on extension
      if (grepl("\\.xlsx?$", file_name, ignore.case = TRUE)) {
        raw <- read_excel(file_path,
                          col_names = input$has_header,
                          na = c("Undetermined", "#VALUE!", "NA", ""))
      } else {
        raw <- read.csv(file_path,
                        header = input$has_header,
                        stringsAsFactors = FALSE,
                        na.strings = c("Undetermined", "#VALUE!", "NA", ""))
      }
      
      # Clean data
      raw <- raw[, colSums(is.na(raw)) < nrow(raw), drop = FALSE]
      raw <- raw[rowSums(is.na(raw)) < ncol(raw), , drop = FALSE]
      
      if (ncol(raw) == 0 || nrow(raw) == 0) {
        showNotification("Uploaded file has no usable data after cleaning.", type = "error")
        rv$file_uploaded <- FALSE
        return()
      }
      
      # Convert to numeric with comprehensive error handling
      conversion_warnings <- 0
      raw_numeric <- as.data.frame(lapply(raw, function(x) {
        original_length <- length(x)
        converted <- tryCatch({
          if (is.numeric(x)) {
            x
          } else {
            as.numeric(as.character(x))
          }
        }, warning = function(w) {
          conversion_warnings <<- conversion_warnings + 1
          rep(NA, original_length)
        }, error = function(e) {
          conversion_warnings <<- conversion_warnings + 1
          rep(NA, original_length)
        })
        return(converted)
      }))
      
      # Notify user if there were conversion issues
      if (conversion_warnings > 0) {
        showNotification(
          paste("Warning: Some non-numeric values were found and converted to NA. Please check your data format."),
          type = "warning",
          duration = 5
        )
      }
      
      rv$raw_data <- raw_numeric
      rv$file_uploaded <- TRUE
      rv$analysis_complete <- FALSE
      
      # Check for outliers
      rv$has_outliers <- check_for_outliers(raw_numeric, reps_num())
      
      showNotification("File uploaded successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      rv$file_uploaded <- FALSE
    })
  })
  
  # Clipboard → textarea (works on https:// or localhost)
  observeEvent(input$paste_clip, {
    session$sendCustomMessage("clipboardRead", list())
  })
  
  observeEvent(input$clip_payload, {
    if (is.character(input$clip_payload) && !startsWith(input$clip_payload, "ERROR:")) {
      updateTextAreaInput(session, "paste_box", value = input$clip_payload)
      showNotification("Pasted from clipboard.", type = "message")
    } else {
      showNotification("Clipboard blocked by browser. Use Ctrl/Cmd+V to paste manually.", type = "warning")
    }
  })
  
  # Parse pasted data and load as raw_data
  observeEvent(input$load_paste, {
    req(input$paste_box)
    tryCatch({
      raw_numeric <- parse_pasted_matrix(
        txt = input$paste_box,
        has_header = isTRUE(input$paste_has_header),
        has_left_labels = isTRUE(input$paste_has_left_labels),
        transpose = isTRUE(input$paste_transpose)
      )
      # Save and mark as uploaded
      rv$raw_data <- raw_numeric
      rv$file_uploaded <- TRUE
      rv$analysis_complete <- FALSE

      # Outlier flag
      rv$has_outliers <- check_for_outliers(raw_numeric, reps_num())
      showNotification("Pasted data loaded!", type = "message")
      updateTabsetPanel(session, "main_tabs", selected = "👀 Data Preview")
    }, error = function(e) {
      showNotification(paste("Paste parsing error:", e$message), type = "error", duration = 7)
    })
  })
  
  # Analysis button observer
  observeEvent(input$analyze, {
    req(rv$raw_data)
    
    withProgress(message = "Analyzing data...", value = 0, {
      
      tryCatch({
        
        incProgress(0.1, detail = "Cleaning data...")
        
        # Step 1: Clean data and remove outliers
        cleaned <- clean_data(rv$raw_data, reps_num(), input$remove_outliers, input$outlier_threshold)
        rv$cleaned_data <- cleaned$data
        rv$outlier_info <- cleaned$outlier_info
        
        # Validate that columns are a multiple of replicate count
        if ((ncol(rv$cleaned_data) %% reps_num()) != 0) {
          showNotification("Column count is not a multiple of the replicate setting.", type = "error")
          return()
        }
        
        incProgress(0.3, detail = "Computing ΔCt values...")
        
        # Step 2: Compute ΔCt with optimal pairing (returns both signs)
        dct_result <- compute_dct(rv$cleaned_data, reps_num())
        rv$dct_data       <- dct_result$dct_pos   # keep if you still want +ΔCt for reference
        rv$neg_dct_data   <- dct_result$dct_neg   # this is what the sheet shows
        
        incProgress(0.5, detail = "Computing expression values...")
        
        # Step 3: Expression = 2^(ref - target) = 2^(-ΔCt) using –ΔCt directly
        neg_mat <- as.matrix(rv$neg_dct_data); storage.mode(neg_mat) <- "double"
        expr_mat <- 2 ^ (neg_mat)
        # force reference block to exactly 1
        expr_mat[, 1:reps_num()] <- 1
        rv$expr_data <- as.data.frame(expr_mat)
        
        incProgress(0.7, detail = "Computing relative expression...")
        
        # Step 4: Relative expression (sheet style) = ratio of replicate means vs control row 1
        rv$rel_means_data <- compute_relative_expression_ratio_of_means(rv$expr_data, reps_num(), control_row = 1)
        
        # (Optional) keep replicate-wise relative table too (what you already had)
        rv$rel_data <- {
          ctrl_row <- 1
          expr_m <- as.matrix(rv$expr_data)
          rel_rep <- sweep(expr_m, 2, expr_m[ctrl_row, ], "/")
          rel_rep[, 1:reps_num()] <- 1
          as.data.frame(rel_rep)
        }
        
        incProgress(0.9, detail = "Preparing output...")
        
        # Step 5: Prepare CSV output
        rv$csv_data <- prepare_csv_output(rv$cleaned_data, rv$neg_dct_data, rv$expr_data, rv$rel_data, reps_num())
        
        # Step 6: Calculate QC metrics
        rv$qc_metrics <- calculate_qc_metrics(rv$cleaned_data, rv$expr_data, rv$rel_data, reps_num())
        
        rv$analysis_complete <- TRUE
        
        incProgress(1, detail = "Analysis complete!")
        
        showNotification("Analysis completed successfully!", type = "message")
        
        # Switch to results tab
        updateTabsetPanel(session, "main_tabs", selected = "📈 Analysis Results")
        
      }, error = function(e) {
        showNotification(paste("Analysis error:", e$message), type = "error")
      })
    })
  })
  
  # Reset button observer
  observeEvent(input$reset, {
    rv$raw_data <- NULL
    rv$cleaned_data <- NULL
    rv$dct_data <- NULL
    rv$expr_data <- NULL
    rv$rel_data <- NULL
    rv$csv_data <- NULL
    rv$analysis_complete <- FALSE
    rv$file_uploaded <- FALSE
    rv$has_outliers <- FALSE
    rv$qc_metrics <- NULL
    rv$outlier_info <- NULL
    
    reset("file")
    updateTabsetPanel(session, "main_tabs", selected = "📖 Instructions")
    
    showNotification("App reset successfully!", type = "message")
  })
  
  # Output functions
  output$file_uploaded <- reactive(rv$file_uploaded)
  output$analysis_complete <- reactive(rv$analysis_complete)
  output$has_outliers <- reactive(rv$has_outliers)
  
  # Unsuspend reactive outputs so UI conditionals update
  outputOptions(output, "file_uploaded", suspendWhenHidden = FALSE)
  outputOptions(output, "analysis_complete", suspendWhenHidden = FALSE)
  outputOptions(output, "has_outliers", suspendWhenHidden = FALSE)
  
  output$file_info <- renderText({
    req(rv$raw_data)
    paste("File loaded successfully!\n",
          "Rows:", nrow(rv$raw_data), "\n",
          "Columns:", ncol(rv$raw_data), "\n",
          "Estimated genes:", floor(ncol(rv$raw_data) / reps_num()))
  })
  
  output$raw_data_table <- renderDT({
    req(rv$raw_data)
    datatable(
      head(rv$raw_data, 10),
      options = list(
        pageLength = 5,
        scrollX = TRUE,
        dom = 't'
      ),
      rownames = FALSE
    )
  })
  
  output$ct_table <- renderDT({
    req(rv$cleaned_data)
    datatable(
      round(rv$cleaned_data, 3),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      caption = "Cleaned Ct Values (outliers removed)"
    )
  })
  
  output$dct_table <- renderDT({
    req(rv$neg_dct_data)
    datatable(
      round(rv$neg_dct_data, 3),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      caption = "−ΔCt (Ct_ref − Ct_target) — sheet style"
    )
  })
  
  output$expr_table <- renderDT({
    req(rv$expr_data)
    datatable(
      round(rv$expr_data, 6),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      caption = "Expression Values (2^ΔCt)"
    )
  })
  
  output$rel_table <- renderDT({
    req(rv$rel_data)
    datatable(
      round(rv$rel_data, 6),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      caption = "Relative Expression vs. Control"
    )
  })

  output$rel_means_table <- renderDT({
    req(rv$rel_means_data)
    datatable(
      round(rv$rel_means_data, 6),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      caption = "Fold change vs control (ratio of replicate means) — sheet style"
    )
  })
  
  output$expr_summary_table <- renderTable({
    req(rv$qc_metrics)
    rv$qc_metrics$expr_summary
  }, rownames = TRUE, digits = 6)
  
  output$rel_summary_table <- renderTable({
    req(rv$qc_metrics)
    rv$qc_metrics$rel_summary
  }, rownames = TRUE, digits = 6)
  
  output$qc_summary <- renderTable({
    req(rv$qc_metrics)
    rv$qc_metrics$qc_summary
  }, rownames = FALSE, digits = 3)
  
  output$outlier_summary <- renderTable({
    req(rv$outlier_info)
    rv$outlier_info
  }, rownames = FALSE, digits = 0)
  

  
  # Download handlers
  output$download_csv_stacked <- downloadHandler(
    filename = function() paste0("qpcr_results_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(rv$csv_data)
      write.table(rv$csv_data,
        file = file,
        sep = ",",
        col.names = FALSE,
        row.names = FALSE,
        quote = TRUE,   # keep section titles safely quoted
        na = ""
      )
    }
  )

  output$download_csv_tidy <- downloadHandler(
    filename = function() paste0("qpcr_tidy_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(rv$cleaned_data, rv$neg_dct_data, rv$expr_data, rv$rel_data)
      reps_now <- as.integer(input$reps)
      tidy <- make_tidy_long(rv$cleaned_data, rv$neg_dct_data, rv$expr_data, rv$rel_data, reps_now)
      write.csv(tidy, file, row.names = FALSE, na = "")
    }
  )
  
  # Helper functions
  check_for_outliers <- function(data, reps) {
    reps <- as.integer(reps)
    if (reps != 4) return(FALSE)
    
    n_genes <- ncol(data) / reps
    has_outliers <- FALSE
    
    for (g in 1:n_genes) {
      start_col <- (g - 1) * reps + 1
      end_col <- g * reps
      
      for (row in seq_len(nrow(data))) {
        block <- tryCatch({
          as.numeric(data[row, start_col:end_col])
        }, warning = function(w) {
          rep(NA, length(start_col:end_col))
        }, error = function(e) {
          rep(NA, length(start_col:end_col))
        })
        
        if (sum(!is.na(block)) >= 3) {
          mean_val <- mean(block, na.rm = TRUE)
          sd_val <- sd(block, na.rm = TRUE)
          if (!is.na(mean_val) && !is.na(sd_val) && any(abs(block - mean_val) > 2 * sd_val, na.rm = TRUE)) {
            has_outliers <- TRUE
            break
          }
        }
      }
      if (has_outliers) break
    }
    
    return(has_outliers)
  }
  
  clean_data <- function(raw_data, reps, remove_outliers, threshold) {
    reps <- as.integer(reps)
    cleaned <- raw_data
    outlier_info <- data.frame(
      Gene = character(),
      Sample = integer(),
      Outliers_Removed = integer(),
      stringsAsFactors = FALSE
    )
    
    if (remove_outliers && reps == 4) {
      n_genes <- ncol(raw_data) / reps
      
      for (g in 1:n_genes) {
        start_col <- (g - 1) * reps + 1
        end_col <- g * reps
        
        for (row in seq_len(nrow(raw_data))) {
          block <- tryCatch({
            as.numeric(raw_data[row, start_col:end_col])
          }, warning = function(w) {
            rep(NA, length(start_col:end_col))
          }, error = function(e) {
            rep(NA, length(start_col:end_col))
          })
          
          if (sum(!is.na(block)) >= 3) {
            mean_val <- mean(block, na.rm = TRUE)
            sd_val <- sd(block, na.rm = TRUE)
            
            if (!is.na(mean_val) && !is.na(sd_val)) {
              outliers <- which(abs(block - mean_val) > threshold * sd_val)
              
              if (length(outliers) > 0 && length(outliers) < length(block)) {
                # Remove the most extreme outlier
                extreme_outlier <- outliers[which.max(abs(block[outliers] - mean_val))]
                cleaned[row, start_col + extreme_outlier - 1] <- NA
                
                outlier_info <- rbind(outlier_info, data.frame(
                  Gene = paste0("Gene", g),
                  Sample = row,
                  Outliers_Removed = 1,
                  stringsAsFactors = FALSE
                ))
              }
            }
          }
        }
      }
    }
    
    return(list(data = cleaned, outlier_info = outlier_info))
  }
  
  compute_dct <- function(cleaned_data, reps) {
    reps <- as.integer(reps)
    n_genes <- ncol(cleaned_data) / reps

    gene_indices <- lapply(0:(n_genes-1), function(i) {
      s <- i * reps + 1; e <- s + reps - 1; s:e
    })

    # minimize SD across replicates by permuting the *target* against ref
    # Returns ΔCt = target - ref (pos_dCt) and -ΔCt = ref - target (neg_dCt)
    compute_pair_min_sd <- function(ref_vals, target_vals) {
      if (all(is.na(ref_vals)) || all(is.na(target_vals))) {
        return(list(pos = rep(NA_real_, length(target_vals)),
                    neg = rep(NA_real_, length(target_vals))))
      }
      # if too sparse, fallback index-wise
      if (sum(!is.na(target_vals)) < 2 || sum(!is.na(ref_vals)) < 2) {
        pos <- target_vals - ref_vals
        pos[is.na(ref_vals) | is.na(target_vals)] <- NA_real_
        return(list(pos = pos, neg = -pos))
      }
      best_pos <- NULL; best_sd <- Inf
      # try permutations; reps ∈ {3,4} so at most 24 perms → cheap
      try({
        perms <- permutations(n = length(target_vals),
                              r = length(target_vals),
                              v = target_vals)
        for (i in seq_len(nrow(perms))) {
          pair <- perms[i, ]
          pos <- pair - ref_vals          # ΔCt = target - ref
          pos[is.na(ref_vals) | is.na(pair)] <- NA_real_
          if (sum(!is.na(pos)) >= 2) {
            s <- sd(pos, na.rm = TRUE)
            if (s < best_sd) { best_sd <- s; best_pos <- pos }
          }
        }
      }, silent = TRUE)
      if (is.null(best_pos)) {
        best_pos <- target_vals - ref_vals
        best_pos[is.na(ref_vals) | is.na(target_vals)] <- NA_real_
      }
      list(pos = best_pos, neg = -best_pos)  # neg = ref - target
    }

    pos_list <- vector("list", nrow(cleaned_data))
    neg_list <- vector("list", nrow(cleaned_data))

    for (r in seq_len(nrow(cleaned_data))) {
      row_vals <- suppressWarnings(as.numeric(cleaned_data[r, ]))
      ref_vals <- row_vals[gene_indices[[1]]]

      # reference gene ΔCt = 0; –ΔCt = 0
      pos_row <- rep(0, reps)
      neg_row <- rep(0, reps)

      for (g in 2:n_genes) {
        target_vals <- row_vals[gene_indices[[g]]]
        pair <- compute_pair_min_sd(ref_vals, target_vals)
        pos_row <- c(pos_row, pair$pos)
        neg_row <- c(neg_row, pair$neg)
      }
      pos_list[[r]] <- pos_row
      neg_list[[r]] <- neg_row
    }

    list(
      dct_pos = as.data.frame(do.call(rbind, pos_list)),   # ΔCt = target - ref
      dct_neg = as.data.frame(do.call(rbind, neg_list))    # –ΔCt = ref - target  (what your sheet shows)
    )
  }
  
  compute_relative_expression_ratio_of_means <- function(expr_data, reps, control_row = 1) {
    reps <- as.integer(reps)
    M <- as.matrix(expr_data); storage.mode(M) <- "double"
    stopifnot(ncol(M) %% reps == 0)
    n_genes <- ncol(M) / reps

    # mean across replicates (per sample × gene block)
    block_mean <- function(m, i1, i2) apply(m[, i1:i2, drop = FALSE], 1, function(x) mean(x, na.rm = TRUE))
    mean_mat <- sapply(seq_len(n_genes), function(g) {
      i1 <- (g-1)*reps + 1; i2 <- g*reps
      block_mean(M, i1, i2)
    })
    ctrl <- mean_mat[control_row, ]
    rel <- sweep(mean_mat, 2, ctrl, "/")    # ratio of means, NOT mean of ratios
    # set reference gene column to 1 explicitly
    rel[, 1] <- 1
    as.data.frame(rel)  # rows = samples, cols = Gene1..GeneN
  }


  
  calculate_qc_metrics <- function(cleaned_data, expr_data, rel_data, reps) {
    reps <- as.integer(reps)
    n_genes <- ncol(cleaned_data) / reps
    
    # Calculate means and SDs
    calc_summary <- function(data, n_reps) {
      means <- list()
      sds <- list()
      cvs <- list()
      
      for (g in 1:n_genes) {
        start_col <- (g - 1) * n_reps + 1
        end_col <- g * n_reps
        
        block_means <- apply(data[, start_col:end_col], 1, function(x) {
          if (all(is.na(x))) NA else mean(x, na.rm = TRUE)
        })
        
        block_sds <- apply(data[, start_col:end_col], 1, function(x) {
          if (all(is.na(x))) NA else sd(x, na.rm = TRUE)
        })
        
        block_cvs <- block_sds / block_means * 100
        
        means[[g]] <- block_means
        sds[[g]] <- block_sds
        cvs[[g]] <- block_cvs
      }
      
      return(list(
        means = do.call(cbind, means),
        sds = do.call(cbind, sds),
        cvs = do.call(cbind, cvs)
      ))
    }
    
    expr_summary <- calc_summary(expr_data, reps)
    rel_summary <- calc_summary(rel_data, reps)
    
    # Create summary tables
    expr_summary_table <- data.frame(
      Gene = paste0("Gene", 1:n_genes),
      Mean = colMeans(expr_summary$means, na.rm = TRUE),
      SD = colMeans(expr_summary$sds, na.rm = TRUE),
      CV_percent = colMeans(expr_summary$cvs, na.rm = TRUE)
    )
    
    rel_summary_table <- data.frame(
      Gene = paste0("Gene", 1:n_genes),
      Mean = colMeans(rel_summary$means, na.rm = TRUE),
      SD = colMeans(rel_summary$sds, na.rm = TRUE),
      CV_percent = colMeans(rel_summary$cvs, na.rm = TRUE)
    )
    
    # QC summary
    qc_summary <- data.frame(
      Metric = c("Total Samples", "Total Genes", "Missing Values (%)", "Mean CV (%)"),
      Value = c(
        nrow(cleaned_data),
        n_genes,
        round(sum(is.na(cleaned_data)) / length(as.matrix(cleaned_data)) * 100, 1),
        round(mean(colMeans(expr_summary$cvs, na.rm = TRUE), na.rm = TRUE), 1)
      )
    )
    
    return(list(
      expr_summary = expr_summary_table,
      rel_summary = rel_summary_table,
      qc_summary = qc_summary,
      cv_data = expr_summary$cvs
    ))
  }
  
  # Helper function to create tidy long format CSV
  make_tidy_long <- function(cleaned_data, dct_data, expr_data, rel_data, reps) {
    reps <- as.integer(reps)
    stopifnot(ncol(cleaned_data) %% reps == 0)

    melt_mat <- function(df, label) {
      m <- as.matrix(df)
      storage.mode(m) <- "double"
      S <- nrow(m); K <- ncol(m)
      gene_idx <- ceiling(seq_len(K) / reps)
      rep_idx  <- ((seq_len(K) - 1) %% reps) + 1
      out <- data.frame(
        Sample = rep(seq_len(S), times = K),
        Gene   = rep(gene_idx, each = S),
        Rep    = rep(rep_idx,  each = S),
        val    = as.vector(m),
        check.names = FALSE
      )
      names(out)[4] <- label
      out
    }

    ct_long   <- melt_mat(cleaned_data, "Ct")
    dct_long  <- melt_mat(dct_data,     "dCt")
    expr_long <- melt_mat(expr_data,    "Expr_2^dCt")
    rel_long  <- melt_mat(rel_data,     "Rel_2^ddCt")

    tidy <- Reduce(function(a,b) merge(a,b, by = c("Sample","Gene","Rep"), all = TRUE),
                   list(ct_long, dct_long, expr_long, rel_long))

    # Human-friendly gene labels (Gene1..GeneN)
    tidy$Gene <- paste0("Gene", tidy$Gene)
    tidy[order(tidy$Sample, tidy$Gene, tidy$Rep), ]
  }

  prepare_csv_output <- function(cleaned_data, dct_data, expr_data, rel_data, reps) {
    reps <- as.integer(reps)
    # This function prepares the comprehensive CSV output
    # Implementation similar to the original but with better formatting
    
    round_data_frame <- function(df, digits = 6) {
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        return(data.frame())
      }
      df[] <- lapply(df, function(x) {
        if (is.numeric(x)) round(x, digits) else x
      })
      return(df)
    }
    
    clean_csv_block <- function(df) {
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        return(matrix("", nrow = 1, ncol = 1))
      }
      mat <- as.matrix(df)
      mat[is.na(mat)] <- "#VALUE!"
      return(apply(mat, c(1,2), as.character))
    }
    
    # Prepare all data blocks
    ct_block <- clean_csv_block(round_data_frame(cleaned_data))
    dct_block <- clean_csv_block(round_data_frame(dct_data))
    expr_block <- clean_csv_block(round_data_frame(expr_data))
    rel_block <- clean_csv_block(round_data_frame(rel_data))
    
    # Calculate summaries
    n_genes <- ncol(cleaned_data) / reps
    
    calc_summary <- function(data, n_reps) {
      means <- list()
      sds <- list()
      
      for (g in 1:n_genes) {
        start_col <- (g - 1) * n_reps + 1
        end_col <- g * n_reps
        
        block_means <- apply(data[, start_col:end_col], 1, function(x) {
          if (all(is.na(x))) NA else mean(x, na.rm = TRUE)
        })
        
        block_sds <- apply(data[, start_col:end_col], 1, function(x) {
          if (all(is.na(x))) NA else sd(x, na.rm = TRUE)
        })
        
        means[[g]] <- block_means
        sds[[g]] <- block_sds
      }
      
      return(list(
        means = do.call(cbind, means),
        sds = do.call(cbind, sds)
      ))
    }
    
    expr_summary <- calc_summary(expr_data, reps)
    rel_summary <- calc_summary(rel_data, reps)
    
    expr_mean_block <- clean_csv_block(round_data_frame(expr_summary$means))
    expr_sd_block <- clean_csv_block(round_data_frame(expr_summary$sds))
    rel_mean_block <- clean_csv_block(round_data_frame(rel_summary$means))
    rel_sd_block <- clean_csv_block(round_data_frame(rel_summary$sds))
    
    # Combine all blocks
    max_cols <- max(
      ncol(ct_block), ncol(dct_block), ncol(expr_block), ncol(rel_block),
      ncol(expr_mean_block), ncol(expr_sd_block), ncol(rel_mean_block), ncol(rel_sd_block),
      na.rm = TRUE
    )
    
    pad_block <- function(block, target_cols) {
      if (ncol(block) < target_cols) {
        extra_cols <- target_cols - ncol(block)
        block <- cbind(block, matrix("", nrow = nrow(block), ncol = extra_cols))
      }
      return(block)
    }
    
    ct_block <- pad_block(ct_block, max_cols)
    dct_block <- pad_block(dct_block, max_cols)
    expr_block <- pad_block(expr_block, max_cols)
    rel_block <- pad_block(rel_block, max_cols)
    expr_mean_block <- pad_block(expr_mean_block, max_cols)
    expr_sd_block <- pad_block(expr_sd_block, max_cols)
    rel_mean_block <- pad_block(rel_mean_block, max_cols)
    rel_sd_block <- pad_block(rel_sd_block, max_cols)
    
    # Create titles
    title_ct <- matrix(c("### Cleaned Ct Data ###", rep("", max_cols - 1)), nrow = 1)
    title_dct <- matrix(c("### −ΔCt (Ct_ref - Ct_target) ###", rep("", max_cols - 1)), nrow = 1)
    title_expr <- matrix(c("### Expression Values (2^(ref - target)) ###", rep("", max_cols - 1)), nrow = 1)
    title_expr_mean <- matrix(c("### Expression Means ###", rep("", max_cols - 1)), nrow = 1)
    title_expr_sd <- matrix(c("### Expression SDs ###", rep("", max_cols - 1)), nrow = 1)
    title_rel <- matrix(c("### Relative Expression vs. Control ###", rep("", max_cols - 1)), nrow = 1)
    title_rel_mean <- matrix(c("### Relative Expression Means ###", rep("", max_cols - 1)), nrow = 1)
    title_rel_sd <- matrix(c("### Relative Expression SDs ###", rep("", max_cols - 1)), nrow = 1)
    
    blank_row <- matrix("", nrow = 1, ncol = max_cols)
    
    final_csv <- rbind(
      title_ct, ct_block, blank_row,
      title_dct, dct_block, blank_row,
      title_expr, expr_block, blank_row,
      title_expr_mean, expr_mean_block, blank_row,
      title_expr_sd, expr_sd_block, blank_row,
      title_rel, rel_block, blank_row,
      title_rel_mean, rel_mean_block, blank_row,
      title_rel_sd, rel_sd_block
    )
    
    return(final_csv)
  }
  

}

shinyApp(ui, server)