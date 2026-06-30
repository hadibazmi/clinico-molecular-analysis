# ---------------------------------------------------------------------------------------
#                                      ROC analysis
# This function performs the entire ROC analysis for a single gene (case and control cols).
# It computes AUC, 95% CI, best threshold (Youden), Sensitivity, Specificity, PPV, NPV,
# LR+, LR-, and also plots and saves the ROC curve in multiple formats.
# ---------------------------------------------------------------------------------------

# ---- Main function ----
analyze_roc_gene = function(data, case_col, control_col, gene_name, 
                            img_width = 7, img_height = 7, dpi_val = 300,
                            # New arguments for precise color control
                            roc_color = "#1976D2", # Default main color (Blue)
                            diag_color = "gray40"   # Default diagonal line color (Gray)
) {
  
  # --- Statistical calculations ---
  # Create a temporary data frame with the two columns and remove rows with any NA values
  temp_df = data.frame(
    Case = data[[case_col]],
    Control = data[[control_col]]
  )
  temp_df = na.omit(temp_df)
  
  # Convert columns to numeric vectors
  case_vals = as.vector(temp_df$Case)
  control_vals = as.vector(temp_df$Control)
  
  # Combine case and control scores and create binary labels (1 for case, 0 for control)
  scores = c(case_vals, control_vals)
  labels = c(rep(1, length(case_vals)), rep(0, length(control_vals)))
  
  # Compute ROC curve using pROC
  roc_obj = roc(labels, scores)
  
  # Calculate AUC and 95% CI (lower and upper bounds)
  auc_val = auc(roc_obj)
  ci_obj = ci.auc(roc_obj)
  ci_lower = ci_obj[1]
  ci_upper = ci_obj[3]
  
  # Find the best threshold using Youden's index
  best_coords = coords(
    roc_obj,
    x = "best",
    ret = c("threshold", "sensitivity", "specificity"),
    best.method = "youden"
  )
  
  # Extract the first row if multiple thresholds have the same index
  best_threshold = as.numeric(best_coords$threshold[1])
  sensitivity_val = as.numeric(best_coords$sensitivity[1])
  specificity_val = as.numeric(best_coords$specificity[1])
  
  # Compute PPV and NPV using approximate formulas based on group sizes
  ppv_val = (sensitivity_val * length(case_vals)) /
    ((sensitivity_val * length(case_vals)) + ((1 - specificity_val) * length(control_vals)))
  npv_val = (specificity_val * length(control_vals)) /
    ((specificity_val * length(control_vals)) + ((1 - sensitivity_val) * length(case_vals)))
  
  # Compute Positive and Negative Likelihood Ratios
  positive_lr = sensitivity_val / (1 - specificity_val)
  negative_lr = (1 - sensitivity_val) / specificity_val
  
  # Create a results data frame
  result_df = data.frame(
    Gene = gene_name,
    AUC = auc_val,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    Best_Threshold = best_threshold,
    Sensitivity = sensitivity_val,
    Specificity = specificity_val,
    PPV = ppv_val,
    NPV = npv_val,
    Positive_LR = positive_lr,
    Negative_LR = negative_lr,
    stringsAsFactors = FALSE
  )
  
  # --- Data preparation for plotting ---
  roc_plot_df = data.frame(
    Specificity = (1 - roc_obj$specificities) * 100,
    Sensitivity = roc_obj$sensitivities * 100
  )
  
  # --- Plotting the ROC curve ---
  roc_plot = ggplot(data = roc_plot_df, aes(x = Specificity, y = Sensitivity)) + 
    # Apply the main user-defined color directly to the line and points
    geom_line(size = 1, color = roc_color) +       
    geom_point(size = 2, color = roc_color) +      
    # Apply the user-defined color to the diagonal line
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = diag_color) +
    labs(
      title = paste("ROC Curve for", gene_name),
      x = "100% - Specificity (%)",
      y = "Sensitivity (%)"
    ) +
    scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
    scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      # Completely remove the legend as requested
      legend.position = "none" 
    )
  print(roc_plot)
  
  # --- Saving files ---
  # Create a folder named after the gene (if it doesn't exist)
  if(!dir.exists(gene_name)) {
    dir.create(gene_name, showWarnings = FALSE)
  }
  
  # Save the ROC plot in multiple formats
  ggsave(filename = paste0(gene_name, "/", gene_name, "_ROC.png"), plot = roc_plot,
         dpi = dpi_val, width = img_width, height = img_height, device = "png")
  ggsave(filename = paste0(gene_name, "/", gene_name, "_ROC.jpeg"), plot = roc_plot,
         dpi = dpi_val, width = img_width, height = img_height, device = "jpeg")
  ggsave(filename = paste0(gene_name, "/", gene_name, "_ROC.tiff"), plot = roc_plot,
         dpi = dpi_val, width = img_width, height = img_height, device = "tiff")
  ggsave(filename = paste0(gene_name, "/", gene_name, "_ROC.pdf"), plot = roc_plot,
         dpi = dpi_val, width = img_width, height = img_height, device = "pdf")
  
  return(result_df)
}



# ============================================
# EXAMPLE USAGE
# ============================================

filepath=read_excel("D:filepath/filename.xlsx")

VariableName= analyze_roc_gene(filepath,
                               case_col ="TumoralRelativeExp" ,
                               control_col ="MarginalRelativeExp" ,
                               gene_name ="P53" )

write_xlsx(VariableName, "VariableName.xlsx")

