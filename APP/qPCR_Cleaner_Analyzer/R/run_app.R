#' Run qPCR Cleaner and Analyzer
#'
#' @description
#' Launches the qPCR Cleaner and Analyzer Shiny application.
#' This function starts the interactive web application for processing
#' and analyzing qPCR data.
#'
#' @param launch.browser Logical. Should the app be launched in the default browser?
#' @param port Integer. Port on which the app should run.
#' @param host Character. Host address on which the app should run.
#'
#' @return
#' Launches the Shiny application in a web browser.
#'
#' @examples
#' \dontrun{
#' # Run the app
#' run_app()
#'
#' # Run on a specific port
#' run_app(port = 8080)
#'
#' # Run without opening browser
#' run_app(launch.browser = FALSE)
#' }
#'
#' @export
run_app <- function(launch.browser = TRUE, port = NULL, host = "127.0.0.1") {
  
  # Check if required packages are installed
  required_packages <- c("shiny", "gtools", "readxl", "DT", "ggplot2", 
                        "dplyr", "tidyr", "shinyjs", "shinyWidgets")
  
  missing_packages <- c()
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "), 
         "\nPlease install them with: install.packages(c('", 
         paste(missing_packages, collapse = "', '"), "'))")
  }
  
  # Get the app directory
  app_dir <- system.file("app", package = "qPCR_Cleaner_Analyzer")
  
  if (app_dir == "") {
    # If not installed as package, use current directory
    app_dir <- getwd()
  }
  
  # Check if app.R exists
  app_file <- file.path(app_dir, "app.R")
  if (!file.exists(app_file)) {
    stop("app.R file not found. Please ensure you're in the correct directory.")
  }
  
  # Launch the app
  shiny::runApp(
    appDir = app_dir,
    launch.browser = launch.browser,
    port = port,
    host = host
  )
}
