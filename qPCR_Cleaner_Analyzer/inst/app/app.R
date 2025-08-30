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
library(ggplot2)
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
          downloadButton("downloadData", 
                        "📥 Download Results", 
                        class = "btn-success",
                        style = "width: 100%;")
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
                )
              )
            ),
            
            # Visualizations
            wellPanel(
              h3("📊 Data Visualizations", style = "color: #2c3e50;"),
              
              fluidRow(
                column(6,
                  plotOutput("ct_boxplot", height = "300px")
                ),
                column(6,
                  plotOutput("expr_boxplot", height = "300px")
                )
              ),
              
              fluidRow(
                column(12,
                  plotOutput("rel_heatmap", height = "400px")
                )
              )
            )
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
              
              hr(),
              
              h5("Replicate Variability (CV%):"),
              plotOutput("cv_plot", height = "300px")
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
      
      # Convert to numeric
      raw_numeric <- as.data.frame(lapply(raw, function(x) {
        if (is.numeric(x)) x else as.numeric(as.character(x))
      }))
      
      rv$raw_data <- raw_numeric
      rv$file_uploaded <- TRUE
      rv$analysis_complete <- FALSE
      
      # Check for outliers
      rv$has_outliers <- check_for_outliers(raw_numeric, input$reps)
      
      showNotification("File uploaded successfully!", type = "success")
      
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      rv$file_uploaded <- FALSE
    })
  })
  
  # Analysis button observer
  observeEvent(input$analyze, {
    req(rv$raw_data, input$file)
    
    withProgress(message = "Analyzing data...", value = 0, {
      
      tryCatch({
        
        incProgress(0.1, detail = "Cleaning data...")
        
        # Step 1: Clean data and remove outliers
        cleaned <- clean_data(rv$raw_data, input$reps, input$remove_outliers, input$outlier_threshold)
        rv$cleaned_data <- cleaned$data
        rv$outlier_info <- cleaned$outlier_info
        
        incProgress(0.3, detail = "Computing ΔCt values...")
        
        # Step 2: Compute ΔCt
        dct_result <- compute_dct(rv$cleaned_data, input$reps)
        rv$dct_data <- dct_result$dct_data
        
        incProgress(0.5, detail = "Computing expression values...")
        
        # Step 3: Compute expression values
        rv$expr_data <- 2 ^ rv$dct_data
        
        incProgress(0.7, detail = "Computing relative expression...")
        
        # Step 4: Compute relative expression
        rv$rel_data <- compute_relative_expression(rv$expr_data, input$reps)
        
        incProgress(0.9, detail = "Preparing output...")
        
        # Step 5: Prepare CSV output
        rv$csv_data <- prepare_csv_output(rv$cleaned_data, rv$dct_data, rv$expr_data, rv$rel_data, input$reps)
        
        # Step 6: Calculate QC metrics
        rv$qc_metrics <- calculate_qc_metrics(rv$cleaned_data, rv$expr_data, rv$rel_data, input$reps)
        
        rv$analysis_complete <- TRUE
        
        incProgress(1, detail = "Analysis complete!")
        
        showNotification("Analysis completed successfully!", type = "success")
        
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
    
    showNotification("App reset successfully!", type = "info")
  })
  
  # Output functions
  output$file_uploaded <- reactive(rv$file_uploaded)
  output$analysis_complete <- reactive(rv$analysis_complete)
  output$has_outliers <- reactive(rv$has_outliers)
  
  output$file_info <- renderText({
    req(rv$raw_data)
    paste("File loaded successfully!\n",
          "Rows:", nrow(rv$raw_data), "\n",
          "Columns:", ncol(rv$raw_data), "\n",
          "Estimated genes:", floor(ncol(rv$raw_data) / input$reps))
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
    req(rv$dct_data)
    datatable(
      round(rv$dct_data, 3),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      caption = "ΔCt Values"
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
  
  # Plots
  output$ct_boxplot <- renderPlot({
    req(rv$cleaned_data)
    plot_ct_boxplot(rv$cleaned_data, input$reps)
  })
  
  output$expr_boxplot <- renderPlot({
    req(rv$expr_data)
    plot_expr_boxplot(rv$expr_data, input$reps)
  })
  
  output$rel_heatmap <- renderPlot({
    req(rv$rel_data)
    plot_rel_heatmap(rv$rel_data, input$reps)
  })
  
  output$cv_plot <- renderPlot({
    req(rv$qc_metrics)
    plot_cv_plot(rv$qc_metrics$cv_data, input$reps)
  })
  
  # Download handler
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("qpcr_analysis_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      write.table(rv$csv_data,
                  file = file,
                  sep = ",",
                  col.names = FALSE,
                  row.names = FALSE,
                  quote = FALSE,
                  na = "")
    }
  )
  
  # Helper functions
  check_for_outliers <- function(data, reps) {
    if (reps != 4) return(FALSE)
    
    n_genes <- ncol(data) / reps
    has_outliers <- FALSE
    
    for (g in 1:n_genes) {
      start_col <- (g - 1) * reps + 1
      end_col <- g * reps
      
      for (row in 1:nrow(data)) {
        block <- as.numeric(data[row, start_col:end_col])
        if (sum(!is.na(block)) >= 3) {
          mean_val <- mean(block, na.rm = TRUE)
          sd_val <- sd(block, na.rm = TRUE)
          if (any(abs(block - mean_val) > 2 * sd_val, na.rm = TRUE)) {
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
        
        for (row in 1:nrow(raw_data)) {
          block <- as.numeric(raw_data[row, start_col:end_col])
          
          if (sum(!is.na(block)) >= 3) {
            mean_val <- mean(block, na.rm = TRUE)
            sd_val <- sd(block, na.rm = TRUE)
            
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
    
    return(list(data = cleaned, outlier_info = outlier_info))
  }
  
  compute_dct <- function(cleaned_data, reps) {
    n_reps <- reps
    n_genes <- ncol(cleaned_data) / n_reps
    
    gene_indices <- lapply(0:(n_genes-1), function(i) {
      start <- i * n_reps + 1
      end <- start + n_reps - 1
      start:end
    })
    
    compute_dCt_optimal <- function(control_vals, target_vals) {
      if (all(is.na(control_vals)) || all(is.na(target_vals))) {
        return(rep(NA, length(target_vals)))
      }
      
      # Simple pairing if not enough data for permutations
      if (sum(!is.na(target_vals)) < 2) {
        dCt <- control_vals - target_vals
        dCt[is.na(control_vals) | is.na(target_vals)] <- NA
        return(dCt)
      }
      
      # Try permutations for optimal pairing
      tryCatch({
        perms <- permutations(n = length(target_vals),
                             r = length(target_vals),
                             v = target_vals)
        
        best_dCt <- NULL
        min_sd <- Inf
        
        for (i in 1:nrow(perms)) {
          pair <- perms[i, ]
          dCt <- control_vals - pair
          dCt[is.na(control_vals) | is.na(pair)] <- NA
          
          if (sum(!is.na(dCt)) >= 2) {
            current_sd <- sd(dCt, na.rm = TRUE)
            if (current_sd < min_sd) {
              min_sd <- current_sd
              best_dCt <- dCt
            }
          }
        }
        
        if (!is.null(best_dCt)) {
          return(best_dCt)
        }
      }, error = function(e) {
        # Fall back to simple pairing
      })
      
      # Fallback: simple pairing
      dCt <- control_vals - target_vals
      dCt[is.na(control_vals) | is.na(target_vals)] <- NA
      return(dCt)
    }
    
    dCt_all <- list()
    
    for (row in 1:nrow(cleaned_data)) {
      row_vals <- as.numeric(cleaned_data[row, ])
      ref_vals <- row_vals[gene_indices[[1]]]
      dCt_row <- c()
      
      for (g in 2:n_genes) {
        target_vals <- row_vals[gene_indices[[g]]]
        dCt_vals <- compute_dCt_optimal(ref_vals, target_vals)
        dCt_row <- c(dCt_row, dCt_vals)
      }
      
      dCt_row <- c(rep(NA, n_reps), dCt_row)
      dCt_all[[row]] <- dCt_row
    }
    
    dCt_df <- as.data.frame(do.call(rbind, dCt_all))
    
    return(list(dct_data = dCt_df))
  }
  
  compute_relative_expression <- function(expr_data, reps) {
    control_expr <- as.numeric(expr_data[1, ])
    
    compute_relative_optimal <- function(control_vals, treat_vals) {
      if (all(is.na(control_vals)) || all(is.na(treat_vals))) {
        return(rep(NA, length(treat_vals)))
      }
      
      # Simple division if not enough data
      if (sum(!is.na(treat_vals)) < 2) {
        rel_vals <- treat_vals / control_vals
        rel_vals[is.na(treat_vals) | is.na(control_vals)] <- NA
        return(rel_vals)
      }
      
      # Try permutations for optimal pairing
      tryCatch({
        perms <- permutations(n = length(control_vals),
                             r = length(control_vals),
                             v = control_vals)
        
        best_rel <- NULL
        min_sd <- Inf
        
        for (i in 1:nrow(perms)) {
          perm_ctrl <- perms[i, ]
          rel_vals <- treat_vals / perm_ctrl
          rel_vals[is.na(treat_vals) | is.na(perm_ctrl)] <- NA
          
          if (sum(!is.na(rel_vals)) >= 2) {
            sd_val <- sd(rel_vals, na.rm = TRUE)
            if (sd_val < min_sd) {
              min_sd <- sd_val
              best_rel <- rel_vals
            }
          }
        }
        
        if (!is.null(best_rel)) {
          return(best_rel)
        }
      }, error = function(e) {
        # Fall back to simple division
      })
      
      # Fallback: simple division
      rel_vals <- treat_vals / control_vals
      rel_vals[is.na(treat_vals) | is.na(control_vals)] <- NA
      return(rel_vals)
    }
    
    rel_all <- list()
    
    for (row in 1:nrow(expr_data)) {
      if (row == 1) {
        rel_all[[row]] <- rep(1, ncol(expr_data))
      } else {
        treat_vals <- as.numeric(expr_data[row, ])
        rel_vals <- compute_relative_optimal(control_expr, treat_vals)
        rel_all[[row]] <- rel_vals
      }
    }
    
    rel_df <- as.data.frame(do.call(rbind, rel_all))
    return(rel_df)
  }
  
  calculate_qc_metrics <- function(cleaned_data, expr_data, rel_data, reps) {
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
  
  prepare_csv_output <- function(cleaned_data, dct_data, expr_data, rel_data, reps) {
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
    title_dct <- matrix(c("### ΔCt Values ###", rep("", max_cols - 1)), nrow = 1)
    title_expr <- matrix(c("### Expression Values (2^ΔCt) ###", rep("", max_cols - 1)), nrow = 1)
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
  
  # Plotting functions
  plot_ct_boxplot <- function(cleaned_data, reps) {
    n_genes <- ncol(cleaned_data) / reps
    
    plot_data <- data.frame()
    for (g in 1:n_genes) {
      start_col <- (g - 1) * reps + 1
      end_col <- g * reps
      
      for (row in 1:nrow(cleaned_data)) {
        values <- as.numeric(cleaned_data[row, start_col:end_col])
        plot_data <- rbind(plot_data, data.frame(
          Gene = paste0("Gene", g),
          Sample = paste0("Sample", row),
          Ct_Value = values
        ))
      }
    }
    
    plot_data <- plot_data[!is.na(plot_data$Ct_Value), ]
    
    ggplot(plot_data, aes(x = Gene, y = Ct_Value, fill = Gene)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(title = "Ct Values Distribution by Gene",
           x = "Gene", y = "Ct Value") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  plot_expr_boxplot <- function(expr_data, reps) {
    n_genes <- ncol(expr_data) / reps
    
    plot_data <- data.frame()
    for (g in 1:n_genes) {
      start_col <- (g - 1) * reps + 1
      end_col <- g * reps
      
      for (row in 1:nrow(expr_data)) {
        values <- as.numeric(expr_data[row, start_col:end_col])
        plot_data <- rbind(plot_data, data.frame(
          Gene = paste0("Gene", g),
          Sample = paste0("Sample", row),
          Expression = values
        ))
      }
    }
    
    plot_data <- plot_data[!is.na(plot_data$Expression), ]
    
    ggplot(plot_data, aes(x = Gene, y = Expression, fill = Gene)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(title = "Expression Values Distribution",
           x = "Gene", y = "Expression (2^ΔCt)") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  plot_rel_heatmap <- function(rel_data, reps) {
    n_genes <- ncol(rel_data) / reps
    
    # Calculate mean relative expression for each gene
    mean_rel <- data.frame()
    for (g in 1:n_genes) {
      start_col <- (g - 1) * reps + 1
      end_col <- g * reps
      
      for (row in 1:nrow(rel_data)) {
        values <- as.numeric(rel_data[row, start_col:end_col])
        mean_val <- mean(values, na.rm = TRUE)
        if (!is.na(mean_val)) {
          mean_rel <- rbind(mean_rel, data.frame(
            Gene = paste0("Gene", g),
            Sample = paste0("Sample", row),
            Relative_Expression = mean_val
          ))
        }
      }
    }
    
    ggplot(mean_rel, aes(x = Gene, y = Sample, fill = Relative_Expression)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                          midpoint = 1, limits = c(0, 2)) +
      labs(title = "Relative Expression Heatmap",
           x = "Gene", y = "Sample", fill = "Relative Expression") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  plot_cv_plot <- function(cv_data, reps) {
    n_genes <- ncol(cv_data)
    
    plot_data <- data.frame()
    for (g in 1:n_genes) {
      values <- cv_data[, g]
      plot_data <- rbind(plot_data, data.frame(
        Gene = paste0("Gene", g),
        CV_percent = values
      ))
    }
    
    plot_data <- plot_data[!is.na(plot_data$CV_percent), ]
    
    ggplot(plot_data, aes(x = Gene, y = CV_percent, fill = Gene)) +
      geom_boxplot(alpha = 0.7) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "red", alpha = 0.7) +
      labs(title = "Coefficient of Variation (CV%) by Gene",
           x = "Gene", y = "CV (%)",
           caption = "Red dashed line indicates 20% CV threshold") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
  }
}

shinyApp(ui, server)