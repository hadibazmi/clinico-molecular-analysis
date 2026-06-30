# ==============================================================================
# qpcr-analyze
# Description: Performs qRT-PCR analysis with Standard Significance Stars (***).
# ==============================================================================

# ---- Main function ----
analyze_qpcr <- function(input_data, case_col, control_col, 
                         save_dir = "qpcr_results", palette = NULL,
                         plot_title = NULL, x_label = NULL, y_label = NULL,
                         case_label = "Case", control_label = "Control") {
  
  # --- STEP 1: LOAD DATA ---
  dataset <- NULL
  if (is.character(input_data)) {
    if (file.exists(input_data)) {
      message(paste("Loading Excel file:", input_data))
      dataset <- read_excel(input_data)
    } else {
      stop(paste("Error: File not found at:", input_data))
    }
  } else if (is.data.frame(input_data)) {
    dataset <- input_data
  } else {
    stop("Error: input_data must be either a file path or a dataframe.")
  }
  
  # --- STEP 2: VALIDATE COLUMNS ---
  if (!case_col %in% names(dataset)) stop(paste("Error: Column '", case_col, "' not found."))
  if (!control_col %in% names(dataset)) stop(paste("Error: Column '", control_col, "' not found."))
  
  # --- STEP 3: PREPARE FOLDER ---
  if (!dir.exists(save_dir)) dir.create(save_dir)
  
  # --- STEP 4: PROCESS DATA ---
  case_dct <- as.numeric(na.omit(dataset[[case_col]]))
  control_dct <- as.numeric(na.omit(dataset[[control_col]]))
  
  # Stats Storage
  stats_results <- data.frame(
    Analysis_Step = character(), Test_Name = character(),
    Statistic_Value = numeric(), P_Value = numeric(),
    Interpretation = character(), stringsAsFactors = FALSE
  )
  
  # --- STEP 5: STATISTICS (T-test on Delta Ct) ---
  # Normality
  shapiro_case <- shapiro.test(case_dct)
  shapiro_control <- shapiro.test(control_dct)
  stats_results[1, ] <- c("Normality Case", "Shapiro-Wilk", round(shapiro_case$statistic,3), shapiro_case$p.value, ifelse(shapiro_case$p.value>0.05,"Normal","Not Normal"))
  stats_results[2, ] <- c("Normality Control", "Shapiro-Wilk", round(shapiro_control$statistic,3), shapiro_control$p.value, ifelse(shapiro_control$p.value>0.05,"Normal","Not Normal"))
  
  # Variance & T-Test
  var_test <- var.test(case_dct, control_dct)
  var_equal <- var_test$p.value > 0.05
  stats_results[3, ] <- c("Variance", "F-test", round(var_test$statistic,3), var_test$p.value, ifelse(var_equal,"Equal","Unequal"))
  
  t_test <- t.test(case_dct, control_dct, var.equal = var_equal)
  t_name <- ifelse(var_equal, "Student's t-test", "Welch's t-test")
  stats_results[4, ] <- c("Hypothesis Test", t_name, round(t_test$statistic,3), t_test$p.value, ifelse(t_test$p.value<0.05,"Significant","Ns"))
  
  # --- STEP 6: FOLD CHANGE ---
  mean_case <- mean(case_dct)
  mean_control <- mean(control_dct)
  
  fc_case_vals <- 2 ^ -(case_dct - mean_control)
  fc_control_vals <- 2 ^ -(control_dct - mean_control)
  
  # Create Dataframe with CUSTOM LABELS
  plot_data_long <- data.frame(
    Group = factor(c(rep(control_label, length(fc_control_vals)), 
                     rep(case_label, length(fc_case_vals))), 
                   levels = c(control_label, case_label)),
    Expression = c(fc_control_vals, fc_case_vals)
  )
  
  # --- STEP 7: PLOTTING ---
  if (is.null(palette)) palette <- c("#999999", "#E69F00") 
  if (is.null(y_label)) y_label <- "Relative Fold Change"
  
  custom_theme <- theme_minimal() +
    theme(
      axis.text = element_text(size = 12, face = "bold", color = "black"),
      axis.title.y = element_text(size = 12, face = "bold", color = "black", margin = margin(r = 10)),
      axis.title.x = element_text(size = 12, face = "bold", color = "black", margin = margin(t = 10)),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 10)),
      plot.subtitle = element_blank()
    )
  
  # --- Determine Significance Stars (CORRECTED LOGIC) ---
  p_val <- t_test$p.value
  sig <- "ns"
  if (p_val < 0.001) {
    sig <- "***"
  } else if (p_val < 0.01) {
    sig <- "**"
  } else if (p_val < 0.05) {
    sig <- "*"
  }
  
  internal_label <- paste0(t_name, "\n(P: ", format.pval(p_val, digits=3), " ", sig, ")")
  
  p <- ggplot(plot_data_long, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black", lwd = 0.6) +
    geom_jitter(width = 0.1, alpha = 0.6, size = 1.8) +
    scale_fill_manual(values = palette) +
    labs(title = plot_title, y = y_label, x = x_label) +
    custom_theme +
    annotate("text", x = 1.5, y = Inf, label = internal_label, 
             vjust = 1.5, size = 3.5, color = "gray40", fontface = "bold")
  
  # --- STEP 8: SAVE ---
  ggsave(filename = file.path(save_dir, "expression_plot_labeled.tiff"), 
         plot = p, width = 6, height = 6, dpi = 300, compression = "lzw")
  
  excel_data <- list("Stats" = stats_results, "PlotData" = plot_data_long, "Raw" = dataset)
  write_xlsx(excel_data, file.path(save_dir, "analysis_results.xlsx"))
  
  message(paste("Done! Results in:", save_dir))
  print(p)
}



# ============================================
# EXAMPLE USAGE
# ============================================

P53=read_excel("filename.xlsx")
analyze_qpcr(input_data = P53,
             case_col = "DeltaCtCase",
             control_col = "DeltaCtControl",
             case_label = "si-P53",
             control_label = "Control",
             y_label = "Expression ratio P53/B2M",
             palette =c("#FF6347","#4682B4") )


