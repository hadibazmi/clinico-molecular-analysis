# ============================================
# Environment Setup Function
# ============================================

init_packages <- function() {
  
  # Define required dependencies
  req_pkgs <- c(
    "readxl", "writexl", "nortest", "ggplot2", 
    "gridExtra", "dplyr", "ggdist", "fs", 
    "openxlsx", "pROC", "pacman", "tidyr"
  )
  
  # Identify installed vs missing packages
  inst_pkgs <- rownames(installed.packages())
  already_inst <- req_pkgs[req_pkgs %in% inst_pkgs]
  missing_pkgs <- req_pkgs[!(req_pkgs %in% inst_pkgs)]
  
  cat("============================================\n")
  cat("          PACKAGE SETUP INITIATED           \n")
  cat("============================================\n\n")
  
  # Report already installed packages
  if (length(already_inst) > 0) {
    cat(">> ALREADY INSTALLED (Skipped):\n")
    cat("   ", paste(already_inst, collapse = ", "), "\n\n")
  }
  
  # Install and report missing packages
  if (length(missing_pkgs) > 0) {
    cat(">> MISSING PACKAGES DETECTED. Installing now:\n")
    cat("   ", paste(missing_pkgs, collapse = ", "), "\n\n")
    
    install.packages(missing_pkgs, repos = "http://cran.us.r-project.org", quiet = TRUE)
    cat(">> Installation of missing packages completed.\n\n")
  } else {
    cat(">> No new installations needed.\n\n")
  }
  
  # Load all required packages silently
  cat(">> Loading packages...\n")
  invisible(suppressPackageStartupMessages(
    lapply(req_pkgs, library, character.only = TRUE)
  ))
  
  # Load core Base R package
  library(grid)
  
  cat(">> Setup successful! Ready for analysis.\n")
  cat("============================================\n")
}

# Setup

init_packages()
