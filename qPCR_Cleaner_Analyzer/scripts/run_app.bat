@echo off
echo ========================================
echo qPCR Cleaner and Analyzer
echo ========================================
echo.
echo Installing required packages and starting the app...
echo.

REM Check if R is installed
where R >nul 2>nul
if %errorlevel% neq 0 (
    echo ERROR: R is not installed or not in PATH
    echo Please install R from https://cran.r-project.org/
    pause
    exit /b 1
)

REM Run the installation script
echo Installing packages...
Rscript install_packages.R

REM Check if installation was successful
if %errorlevel% neq 0 (
    echo.
    echo ERROR: Package installation failed
    echo Please check the error messages above
    pause
    exit /b 1
)

echo.
echo Starting the qPCR Cleaner and Analyzer...
echo The app will open in your default web browser
echo Press Ctrl+C to stop the app
echo.

REM Start the Shiny app
Rscript -e "shiny::runApp(launch.browser = TRUE)"

pause
