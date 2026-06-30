# ============================================
# Demographic Data Extraction Tool (DataFrame Input)
# ============================================

# ---- Main function ----
extract_demographics <- function(dataset, selected_columns, output_name = "Demographic_Table.xlsx") {
  
  # 1. Validate Input
  if (!is.data.frame(dataset)) {
    stop("Error: The input 'dataset' must be a data frame (table).")
  }
  
  data <- dataset
  
  # Initialize list to store results
  table_list <- list()
  counter <- 1
  
  # 2. Loop through user-selected columns
  for (col in selected_columns) {
    
    # Check if column exists
    if (!col %in% names(data)) {
      warning(sprintf("Column '%s' NOT found in data. Skipped.", col))
      next
    }
    
    # Get the data vector
    vec <- data[[col]]
    
    # 3. Determine Data Type (Smart Detection)
    if (is.numeric(vec)) {
      # --- Numeric Analysis (Mean, SD, Range) ---
      vec_clean <- vec[!is.na(vec)]
      
      # Check if data is empty after removing NAs
      if (length(vec_clean) == 0) {
        mean_val <- NA; sd_val <- NA; min_val <- NA; max_val <- NA
      } else {
        mean_val <- mean(vec_clean)
        sd_val   <- sd(vec_clean)
        min_val  <- min(vec_clean)
        max_val  <- max(vec_clean)
      }
      
      row_mean <- data.frame(
        Variable = col,
        Category = "Mean ± SD",
        Count    = sprintf("%.2f ± %.2f", mean_val, sd_val),
        Percent  = "-"
      )
      
      row_range <- data.frame(
        Variable = "",
        Category = "Range (Min - Max)",
        Count    = sprintf("%.2f - %.2f", min_val, max_val),
        Percent  = "-"
      )
      
      table_list[[counter]] <- rbind(row_mean, row_range)
      
    } else {
      # --- Categorical Analysis (Frequency, %) ---
      vec_clean <- as.character(vec)
      vec_clean <- vec_clean[!is.na(vec_clean)]
      total_n   <- length(vec_clean)
      
      tbl <- table(vec_clean)
      df_res <- as.data.frame(tbl)
      colnames(df_res) <- c("Category", "Freq")
      
      df_res$Pct <- (df_res$Freq / total_n) * 100
      
      formatted_rows <- data.frame(
        Variable = c(col, rep("", nrow(df_res) - 1)), 
        Category = df_res$Category,
        Count    = as.character(df_res$Freq),
        Percent  = sprintf("%.1f%%", df_res$Pct)
      )
      
      table_list[[counter]] <- formatted_rows
    }
    
    counter <- counter + 1
  }
  
  # 4. Combine and Save
  if (length(table_list) > 0) {
    final_table <- do.call(rbind, table_list)
    
    # Add Header Row
    final_table <- rbind(
      data.frame(Variable = "Variable", Category = "Sub-group", Count = "N / Value", Percent = "%"),
      final_table
    )
    
    write_xlsx(final_table, output_name)
    message(sprintf(">> Done! Demographic data extracted to: %s", normalizePath(output_name)))
    return(final_table)
    
  } else {
    stop("No valid columns were processed.")
  }
}


# ============================================
# EXAMPLE USAGE
# ============================================

my_file =read_excel("filename.xlsx")

my_vars <- c("Age", "Histology", "Differentiation", "Tumor-size",
             "Lymph-node", "Depth-of-cervical-invasion","Stage","HPV") #Column name in your filename.xlsx

extract_demographics(my_file, my_vars)

# Note: Excel columns representing categorical data but stored as numerics 
# must be explicitly coerced to factors in R prior to downstream analysis.

# Convert categorical variables imported as numerics from Excel into factors for proper extraction.
data$TreatmentGroup <- as.factor(data$TreatmentGroup)














